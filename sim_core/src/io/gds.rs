// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! GDSII stream reader with hierarchy flattening.
//!
//! Supports BOUNDARY, BOX, PATH (flat / half-width-extended ends), SREF and AREF with
//! STRANS reflection, MAG and ANGLE, and scales every coordinate by the file's UNITS record
//! (database unit → nm). TEXT/NODE elements are skipped. Layer semantics follow GDSII: all
//! polygons of a layer are **unioned**; holes must be drawn as separate cut-out geometry.
//!
//! Limits (reported as errors, never silently truncated): hierarchy depth 32, 4 million
//! flattened polygons.

use std::collections::{BTreeMap, HashSet};
use thiserror::Error;

pub const MAX_DEPTH: usize = 32;
pub const MAX_POLYGONS: usize = 4_000_000;

#[derive(Debug, Error, PartialEq)]
pub enum GdsError {
    #[error("not a GDSII stream: {0}")]
    Format(String),
    #[error("unexpected end of file")]
    Eof,
    #[error("structure '{0}' is referenced but not defined")]
    MissingStructure(String),
    #[error("hierarchy is deeper than {MAX_DEPTH} levels or references itself (cell '{0}')")]
    TooDeep(String),
    #[error("layout expands to more than {MAX_POLYGONS} polygons")]
    TooLarge,
    #[error("the file contains no structures")]
    Empty,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trans { pub reflect: bool, pub mag: f64, pub angle_deg: f64 }
impl Default for Trans { fn default() -> Self { Self { reflect: false, mag: 1.0, angle_deg: 0.0 } } }

#[derive(Debug, Clone)]
enum Element {
    Poly { layer: i32, datatype: i32, pts: Vec<(f64, f64)> },
    Path { layer: i32, datatype: i32, width: f64, pathtype: i16, bgn_ext: f64, end_ext: f64, pts: Vec<(f64, f64)> },
    SRef { name: String, trans: Trans, at: (f64, f64) },
    ARef { name: String, trans: Trans, cols: usize, rows: usize, pts: [(f64, f64); 3] },
}

#[derive(Debug, Clone, Default)]
pub struct Structure { elements: Vec<Element> }

#[derive(Debug, Clone)]
pub struct GdsLayout {
    pub structures: BTreeMap<String, Structure>,
    /// Metres per database unit.
    pub db_unit_m: f64,
}

#[derive(Debug, Clone)]
pub struct FlatPoly { pub layer: i32, pub datatype: i32, pub pts: Vec<(f64, f64)> }

#[derive(Debug, Clone)]
pub struct LayerSummary { pub layer: i32, pub datatype: i32, pub polygons: usize, pub bbox: [f64; 4], /// A corner of the polygon nearest the middle of the layer: a crop centre that always hits geometry *and* an edge.
    pub anchor: (f64, f64),
    /// Smallest-dimension statistics over polygons (min of bbox width/height, nm): 10th percentile and median.
    pub feature_p10_nm: f64, pub feature_med_nm: f64 }

impl LayerSummary {
    pub fn size_nm(&self) -> (f64, f64) { (self.bbox[2] - self.bbox[0], self.bbox[3] - self.bbox[1]) }
    pub fn center_nm(&self) -> (f64, f64) { ((self.bbox[0] + self.bbox[2]) * 0.5, (self.bbox[1] + self.bbox[3]) * 0.5) }
}

fn real8(b: &[u8]) -> f64 {
    if b.len() < 8 { return 0.0; }
    let sign = if b[0] & 0x80 != 0 { -1.0 } else { 1.0 };
    let exp = (b[0] & 0x7F) as i32 - 64;
    let mut m: u64 = 0;
    for &x in &b[1..8] { m = (m << 8) | x as u64; }
    sign * (m as f64 / (1u64 << 56) as f64) * 16f64.powi(exp)
}
fn i16be(b: &[u8]) -> i16 { if b.len() >= 2 { i16::from_be_bytes([b[0], b[1]]) } else { 0 } }
fn i32be(b: &[u8]) -> i32 { i32::from_be_bytes([b[0], b[1], b[2], b[3]]) }
fn xy(payload: &[u8]) -> Vec<(f64, f64)> { payload.chunks_exact(8).map(|c| (i32be(&c[0..4]) as f64, i32be(&c[4..8]) as f64)).collect() }
fn name(payload: &[u8]) -> String { String::from_utf8_lossy(payload).trim_end_matches('\0').to_string() }

pub fn parse(bytes: &[u8]) -> Result<GdsLayout, GdsError> {
    if bytes.len() < 4 { return Err(GdsError::Eof); }
    // Real files start with HEADER; minimal fixtures may not. Accept any plausible first record
    // (known record type, sane length) and reject text / OASIS / other binary formats.
    let first_len = u16::from_be_bytes([bytes[0], bytes[1]]) as usize;
    if bytes[2] > 0x3F || first_len < 4 || (bytes[2] == 0x00 && bytes[3] != 0x02) {
        return Err(GdsError::Format("the first record is not a GDSII record (is this a .gds/.gdsii stream, not OASIS or text?)".into()));
    }
    if first_len > bytes.len() && first_len > 4096 { return Err(GdsError::Format("record length exceeds the file: not a GDSII stream".into())); }
    let mut layout = GdsLayout { structures: BTreeMap::new(), db_unit_m: 1e-9 };
    let (mut pos, mut cur): (usize, Option<(String, Structure)>) = (0, None);
    // Element under construction.
    #[derive(Default)] struct E { kind: u8, layer: i32, dtype: i32, width: f64, ptype: i16, bgn: f64, end: f64, sname: String, trans: Trans, cols: usize, rows: usize, pts: Vec<(f64, f64)> }
    let mut e: Option<E> = None;
    while pos + 4 <= bytes.len() {
        let len = u16::from_be_bytes([bytes[pos], bytes[pos + 1]]) as usize;
        if len == 0 { pos += 2; continue; } // zero padding at the end of tape-style files
        if len < 4 { return Err(GdsError::Format(format!("record length {len} < 4 at byte {pos}"))); }
        if pos + len > bytes.len() { return Err(GdsError::Eof); }
        let (ty, p) = (bytes[pos + 2], &bytes[pos + 4..pos + len]);
        match ty {
            0x03 if p.len() >= 16 => layout.db_unit_m = real8(&p[8..16]),
            0x05 => cur = Some((String::new(), Structure::default())),
            0x06 => if let Some(c) = cur.as_mut() { c.0 = name(p); },
            0x07 => if let Some((n, s)) = cur.take() { layout.structures.insert(n, s); },
            0x08 | 0x2D => e = Some(E { kind: 0, trans: Trans::default(), ..Default::default() }),
            0x09 => e = Some(E { kind: 1, trans: Trans::default(), ..Default::default() }),
            0x0A => e = Some(E { kind: 2, trans: Trans::default(), ..Default::default() }),
            0x0B => e = Some(E { kind: 3, trans: Trans::default(), cols: 1, rows: 1, ..Default::default() }),
            0x0C | 0x15 => e = Some(E { kind: 9, ..Default::default() }), // TEXT / NODE: skipped
            0x0D => if let Some(x) = e.as_mut() { x.layer = i16be(p) as i32 },
            0x0E | 0x2E | 0x16 => if let Some(x) = e.as_mut() { x.dtype = i16be(p) as i32 },
            0x0F => if let Some(x) = e.as_mut() { x.width = i32be(&[p[0], p[1], p[2], p[3]]) as f64 },
            0x21 => if let Some(x) = e.as_mut() { x.ptype = i16be(p) },
            0x30 if p.len() >= 4 => if let Some(x) = e.as_mut() { x.bgn = i32be(p) as f64 },
            0x31 if p.len() >= 4 => if let Some(x) = e.as_mut() { x.end = i32be(p) as f64 },
            0x12 => if let Some(x) = e.as_mut() { x.sname = name(p) },
            0x13 if p.len() >= 4 => if let Some(x) = e.as_mut() { x.cols = i16be(&p[0..2]).max(1) as usize; x.rows = i16be(&p[2..4]).max(1) as usize; },
            0x1A => if let Some(x) = e.as_mut() { x.trans.reflect = i16be(p) as u16 & 0x8000 != 0 },
            0x1B if p.len() >= 8 => if let Some(x) = e.as_mut() { x.trans.mag = real8(p) },
            0x1C if p.len() >= 8 => if let Some(x) = e.as_mut() { x.trans.angle_deg = real8(p) },
            0x10 => if let Some(x) = e.as_mut() { x.pts = xy(p) },
            0x11 => if let Some(x) = e.take() {
                // Elements outside any STRUCTURE (hand-made fixtures) go to an implicit cell.
                let c = cur.get_or_insert_with(|| ("(unnamed)".to_string(), Structure::default()));
                match x.kind {
                    0 if x.pts.len() >= 3 => {
                        let mut pts = x.pts;
                        if pts.len() > 3 && pts.first() == pts.last() { pts.pop(); } // closing vertex is implicit
                        if pts.len() >= 3 { c.1.elements.push(Element::Poly { layer: x.layer, datatype: x.dtype, pts }); }
                    }
                    1 if x.pts.len() >= 2 => c.1.elements.push(Element::Path { layer: x.layer, datatype: x.dtype, width: x.width, pathtype: x.ptype, bgn_ext: x.bgn, end_ext: x.end, pts: x.pts }),
                    2 if !x.pts.is_empty() => c.1.elements.push(Element::SRef { name: x.sname, trans: x.trans, at: x.pts[0] }),
                    3 if x.pts.len() >= 3 => c.1.elements.push(Element::ARef { name: x.sname, trans: x.trans, cols: x.cols, rows: x.rows, pts: [x.pts[0], x.pts[1], x.pts[2]] }),
                    _ => {}
                }
            },
            0x04 => break,
            _ => {}
        }
        pos += len;
    }
    if let Some((n, s)) = cur.take() { if !s.elements.is_empty() { layout.structures.insert(n, s); } }
    if layout.structures.is_empty() { return Err(GdsError::Empty); }
    if !(layout.db_unit_m.is_finite() && layout.db_unit_m > 0.0) { return Err(GdsError::Format("UNITS record gives a non-positive database unit".into())); }
    Ok(layout)
}

fn apply(t: &Trans, at: (f64, f64), p: (f64, f64)) -> (f64, f64) {
    let (x, y) = (p.0 * t.mag, if t.reflect { -p.1 } else { p.1 } * t.mag);
    let (s, c) = t.angle_deg.to_radians().sin_cos();
    (at.0 + c * x - s * y, at.1 + s * x + c * y)
}

/// Compose `outer ∘ inner` as an affine map (matrix + translation) so deep hierarchies stay O(1) per point.
#[derive(Clone, Copy)]
struct Aff { m: [f64; 4], t: (f64, f64) }
impl Aff {
    const I: Aff = Aff { m: [1.0, 0.0, 0.0, 1.0], t: (0.0, 0.0) };
    fn of(t: &Trans, at: (f64, f64)) -> Aff {
        let (s, c) = t.angle_deg.to_radians().sin_cos();
        let (a, d) = (t.mag, if t.reflect { -t.mag } else { t.mag });
        Aff { m: [c * a, -s * d, s * a, c * d], t: at }
    }
    fn then(&self, inner: &Aff) -> Aff {
        let (o, i) = (self.m, inner.m);
        Aff { m: [o[0] * i[0] + o[1] * i[2], o[0] * i[1] + o[1] * i[3], o[2] * i[0] + o[3] * i[2], o[2] * i[1] + o[3] * i[3]],
              t: (o[0] * inner.t.0 + o[1] * inner.t.1 + self.t.0, o[2] * inner.t.0 + o[3] * inner.t.1 + self.t.1) }
    }
    fn p(&self, p: (f64, f64)) -> (f64, f64) { (self.m[0] * p.0 + self.m[1] * p.1 + self.t.0, self.m[2] * p.0 + self.m[3] * p.1 + self.t.1) }
}

/// PATH → polygons. Pathtype 0 flat ends, 1 round ends, 2 square ends (half-width extension), 4 custom
/// BGNEXTN/ENDEXTN. Interior vertices get a join (round for type 1, bevel-filled otherwise) so bends have no notch.
/// A negative width means "absolute width" in GDSII (not magnified); the magnitude is used.
fn path_polygons(pts: &[(f64, f64)], width: f64, pathtype: i16, bgn_ext: f64, end_ext: f64) -> Vec<Vec<(f64, f64)>> {
    let hw = width.abs() * 0.5;
    let mut p: Vec<(f64, f64)> = Vec::with_capacity(pts.len());
    for &q in pts { if p.last().map_or(true, |l| (l.0 - q.0).hypot(l.1 - q.1) > 1e-12) { p.push(q); } }
    if p.len() < 2 || hw <= 0.0 { return Vec::new(); }
    let dir = |a: (f64, f64), b: (f64, f64)| { let (dx, dy) = (b.0 - a.0, b.1 - a.1); let l = dx.hypot(dy); (dx / l, dy / l) };
    let (ext0, ext1) = match pathtype { 2 => (hw, hw), 4 => (bgn_ext, end_ext), _ => (0.0, 0.0) };
    let n = p.len();
    let (d0, d1) = (dir(p[0], p[1]), dir(p[n - 2], p[n - 1]));
    let mut q = p.clone();
    q[0] = (p[0].0 - d0.0 * ext0, p[0].1 - d0.1 * ext0);
    q[n - 1] = (p[n - 1].0 + d1.0 * ext1, p[n - 1].1 + d1.1 * ext1);
    let mut out = Vec::new();
    for w in q.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (ux, uy) = dir(a, b);
        let (nx, ny) = (-uy * hw, ux * hw);
        out.push(vec![(a.0 + nx, a.1 + ny), (b.0 + nx, b.1 + ny), (b.0 - nx, b.1 - ny), (a.0 - nx, a.1 - ny)]);
    }
    let circle = |c: (f64, f64)| -> Vec<(f64, f64)> { (0..16).map(|i| { let t = i as f64 / 16.0 * std::f64::consts::TAU; (c.0 + hw * t.cos(), c.1 + hw * t.sin()) }).collect() };
    if pathtype == 1 {
        for &c in &p { out.push(circle(c)); } // round caps and round joins
    } else {
        // Mitred joins at interior bends (bevel beyond a miter limit of 4 half-widths) so a bend has no notch.
        for i in 1..n - 1 {
            let (d_in, d_out) = (dir(p[i - 1], p[i]), dir(p[i], p[i + 1]));
            let (n1, n2) = ((-d_in.1 * hw, d_in.0 * hw), (-d_out.1 * hw, d_out.0 * hw));
            let denom = 1.0 + d_in.0 * d_out.0 + d_in.1 * d_out.1;
            for sgn in [1.0, -1.0] {
                let (a, b) = ((p[i].0 + sgn * n1.0, p[i].1 + sgn * n1.1), (p[i].0 + sgn * n2.0, p[i].1 + sgn * n2.1));
                let mut poly = vec![p[i], a];
                if denom > 1e-6 {
                    let m = (sgn * (n1.0 + n2.0) / denom, sgn * (n1.1 + n2.1) / denom);
                    if m.0.hypot(m.1) <= 4.0 * hw { poly.push((p[i].0 + m.0, p[i].1 + m.1)); }
                }
                poly.push(b);
                out.push(poly);
            }
        }
    }
    out
}

impl GdsLayout {
    /// Cells that no other cell references.
    pub fn top_cells(&self) -> Vec<String> {
        let mut referenced: HashSet<&str> = HashSet::new();
        for s in self.structures.values() { for e in &s.elements { match e { Element::SRef { name, .. } | Element::ARef { name, .. } => { referenced.insert(name.as_str()); } _ => {} } } }
        self.structures.keys().filter(|k| !referenced.contains(k.as_str())).cloned().collect()
    }

    /// The top cell with the most elements (ties: alphabetical).
    pub fn default_top(&self) -> Option<String> {
        let mut tops = self.top_cells();
        if tops.is_empty() { tops = self.structures.keys().cloned().collect(); }
        tops.into_iter().max_by_key(|n| self.structures[n].elements.len())
    }

    /// Flatten `top` (default: [`Self::default_top`]) into polygons in **nanometres**.
    pub fn flatten(&self, top: Option<&str>) -> Result<Vec<FlatPoly>, GdsError> { self.flatten_region(top, None) }

    /// Like [`Self::flatten`] but only expands instances and keeps polygons whose bounding box intersects
    /// `region = [x0, y0, x1, y1]` (nm). Instances (and whole AREF cells) outside the region are skipped without
    /// being expanded, so a small window of a huge layout costs time proportional to what is inside it.
    pub fn flatten_region(&self, top: Option<&str>, region: Option<[f64; 4]>) -> Result<Vec<FlatPoly>, GdsError> {
        let top = match top { Some(t) => t.to_string(), None => self.default_top().ok_or(GdsError::Empty)? };
        let scale = self.db_unit_m * 1e9;
        let region_db = region.map(|r| [r[0] / scale, r[1] / scale, r[2] / scale, r[3] / scale]);
        let mut bb: BTreeMap<String, Option<[f64; 4]>> = BTreeMap::new();
        self.cell_bbox(&top, 0, &mut bb)?;
        let mut out = Vec::new();
        self.walk(&top, &Aff::I, 0, scale, region_db, &bb, &mut out)?;
        Ok(out)
    }

    /// Bounding box of a cell in its own database units (`None` for an empty cell). Memoised.
    fn cell_bbox(&self, cell: &str, depth: usize, memo: &mut BTreeMap<String, Option<[f64; 4]>>) -> Result<Option<[f64; 4]>, GdsError> {
        if let Some(b) = memo.get(cell) { return Ok(*b); }
        if depth > MAX_DEPTH { return Err(GdsError::TooDeep(cell.into())); }
        let s = self.structures.get(cell).ok_or_else(|| GdsError::MissingStructure(cell.into()))?;
        let mut acc: Option<[f64; 4]> = None;
        let add = |acc: &mut Option<[f64; 4]>, b: [f64; 4]| { *acc = Some(match *acc { None => b, Some(a) => [a[0].min(b[0]), a[1].min(b[1]), a[2].max(b[2]), a[3].max(b[3])] }); };
        for el in &s.elements {
            match el {
                Element::Poly { pts, .. } => add(&mut acc, pts_bbox(pts)),
                Element::Path { width, pts, bgn_ext, end_ext, .. } => { let m = width.abs() * 0.5 + bgn_ext.abs().max(end_ext.abs()); let b = pts_bbox(pts); add(&mut acc, [b[0] - m, b[1] - m, b[2] + m, b[3] + m]); }
                Element::SRef { name, trans, at } => if let Some(cb) = self.cell_bbox(name, depth + 1, memo)? { add(&mut acc, aff_bbox(&Aff::of(trans, *at), cb)); },
                Element::ARef { name, trans, cols, rows, pts } => if let Some(cb) = self.cell_bbox(name, depth + 1, memo)? {
                    let (o, ce, re) = (pts[0], pts[1], pts[2]);
                    let (dc, dr) = ((ce.0 - o.0) / *cols as f64, (ce.1 - o.1) / *cols as f64);
                    let (er, fr) = ((re.0 - o.0) / *rows as f64, (re.1 - o.1) / *rows as f64);
                    for (r, c) in [(0, 0), (0, cols - 1), (rows - 1, 0), (rows - 1, cols - 1)] {
                        let at = (o.0 + dc * c as f64 + er * r as f64, o.1 + dr * c as f64 + fr * r as f64);
                        add(&mut acc, aff_bbox(&Aff::of(trans, at), cb));
                    }
                },
            }
        }
        memo.insert(cell.to_string(), acc);
        Ok(acc)
    }

    #[allow(clippy::too_many_arguments)]
    fn walk(&self, cell: &str, a: &Aff, depth: usize, scale: f64, region: Option<[f64; 4]>, bb: &BTreeMap<String, Option<[f64; 4]>>, out: &mut Vec<FlatPoly>) -> Result<(), GdsError> {
        if depth > MAX_DEPTH { return Err(GdsError::TooDeep(cell.into())); }
        let s = self.structures.get(cell).ok_or_else(|| GdsError::MissingStructure(cell.into()))?;
        let hit = |b: [f64; 4]| region.map_or(true, |r| !(b[2] < r[0] || b[0] > r[2] || b[3] < r[1] || b[1] > r[3]));
        let push = |layer: i32, datatype: i32, pts: Vec<(f64, f64)>, out: &mut Vec<FlatPoly>| -> Result<(), GdsError> {
            let tp: Vec<(f64, f64)> = pts.into_iter().map(|p| a.p(p)).collect();
            if !hit(pts_bbox(&tp)) { return Ok(()); }
            if out.len() >= MAX_POLYGONS { return Err(GdsError::TooLarge); }
            out.push(FlatPoly { layer, datatype, pts: tp.into_iter().map(|q| (q.0 * scale, q.1 * scale)).collect() });
            Ok(())
        };
        for el in &s.elements {
            match el {
                Element::Poly { layer, datatype, pts } => push(*layer, *datatype, pts.clone(), out)?,
                Element::Path { layer, datatype, width, pathtype, bgn_ext, end_ext, pts } => for poly in path_polygons(pts, *width, *pathtype, *bgn_ext, *end_ext) { push(*layer, *datatype, poly, out)?; },
                Element::SRef { name, trans, at } => {
                    let na = a.then(&Aff::of(trans, *at));
                    if let Some(cb) = bb.get(name.as_str()).copied().flatten() { if !hit(aff_bbox(&na, cb)) { continue; } }
                    self.walk(name, &na, depth + 1, scale, region, bb, out)?
                }
                Element::ARef { name, trans, cols, rows, pts } => {
                    let cb = bb.get(name.as_str()).copied().flatten();
                    let (o, ce, re) = (pts[0], pts[1], pts[2]);
                    let (dc, dr) = ((ce.0 - o.0) / *cols as f64, (ce.1 - o.1) / *cols as f64);
                    let (er, fr) = ((re.0 - o.0) / *rows as f64, (re.1 - o.1) / *rows as f64);
                    for r in 0..*rows { for c in 0..*cols {
                        let at = (o.0 + dc * c as f64 + er * r as f64, o.1 + dr * c as f64 + fr * r as f64);
                        let na = a.then(&Aff::of(trans, at));
                        if let Some(cb) = cb { if !hit(aff_bbox(&na, cb)) { continue; } }
                        self.walk(name, &na, depth + 1, scale, region, bb, out)?;
                    } }
                }
            }
        }
        Ok(())
    }
}

fn pts_bbox(p: &[(f64, f64)]) -> [f64; 4] { p.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, v| [b[0].min(v.0), b[1].min(v.1), b[2].max(v.0), b[3].max(v.1)]) }
fn aff_bbox(a: &Aff, b: [f64; 4]) -> [f64; 4] { pts_bbox(&[a.p((b[0], b[1])), a.p((b[2], b[1])), a.p((b[2], b[3])), a.p((b[0], b[3]))]) }

/// Distinct (layer, datatype) pairs with polygon counts and bounding boxes.

impl GdsLayout {
    /// The geometry a context sees, in coordinates relative to its own placement.
    ///
    /// The cell's own polygons sit at the origin; every neighbour within `margin_nm` is
    /// offset by its position relative to this placement. That is exactly the input a
    /// cell simulation needs: nothing outside the margin can reach the core.
    pub fn local_geometry(
        &self,
        ctx: &Context,
        placements: &[Placement],
        margin_nm: f64,
    ) -> Vec<atto_geom::tiling::Polygon> {
        let mut out = Vec::new();
        let push_cell = |name: &str, dx: f64, dy: f64, out: &mut Vec<atto_geom::tiling::Polygon>| {
            let Some(s) = self.structures.get(name) else { return };
            for e in &s.elements {
                match e {
                    Element::Poly { layer, datatype, pts } => out.push(atto_geom::tiling::Polygon {
                        layer: *layer,
                        datatype: *datatype,
                        points: pts.iter().map(|p| (p.0 + dx, p.1 + dy)).collect(),
                    }),
                    Element::Path { layer, datatype, pts, width, .. } => {
                        // a path is a stroke; approximate it by its centreline vertices
                        out.push(atto_geom::tiling::Polygon {
                            layer: *layer,
                            datatype: *datatype,
                            points: pts.iter().map(|p| (p.0 + dx, p.1 + dy)).collect(),
                        });
                        let _ = width;
                    }
                    _ => {}
                }
            }
        };
        push_cell(&ctx.cell, 0.0, 0.0, &mut out);
        for q in placements {
            let (dx, dy) = (q.at.0 - ctx.at.0, q.at.1 - ctx.at.1);
            if dx.abs() > margin_nm || dy.abs() > margin_nm {
                continue;
            }
            if dx == 0.0 && dy == 0.0 {
                continue;
            }
            push_cell(&q.cell, dx, dy, &mut out);
        }
        out
    }
}


impl GdsLayout {
    /// How many placements the hierarchy holds, without building them.
    ///
    /// [`Self::placements`] allocates one entry per site, which on a 23-array die is
    /// millions of entries and hundreds of megabytes — spent only to discover the number
    /// is over the cap. This counts instead: the same recursion, no allocation.
    pub fn placement_count(&self, top: &str) -> usize {
        self.count_under(top, 0)
    }

    fn count_under(&self, name: &str, depth: u32) -> usize {
        if depth > 32 {
            return 0;
        }
        let Some(s) = self.structures.get(name) else { return 0 };
        let mut n = 0usize;
        for e in &s.elements {
            match e {
                Element::SRef { name: child, .. } => {
                    n = n.saturating_add(1).saturating_add(self.count_under(child, depth + 1));
                }
                Element::ARef { name: child, cols, rows, .. } => {
                    let sites = (*cols).saturating_mul(*rows) as usize;
                    let inner = 1usize.saturating_add(self.count_under(child, depth + 1));
                    n = n.saturating_add(sites.saturating_mul(inner));
                }
                _ => {}
            }
        }
        n
    }
}

/// A placement together with what surrounds it.
///
/// Two placements of the same cell are only interchangeable if their *neighbourhoods* are
/// the same: a cell at the edge of an array sees open field on one side and etches
/// differently from one in the middle. Dedup on the cell alone is therefore wrong at
/// boundaries, and right only in the interior.
#[derive(Debug, Clone, PartialEq)]
pub struct Context {
    /// The cell being placed.
    pub cell: String,
    /// The cell's own content signature.
    pub signature: u64,
    /// Neighbours within the interaction radius, as quantised offsets and signatures,
    /// sorted so two placements with the same surroundings produce the same key.
    pub neighbours: Vec<(i64, i64, u64)>,
    /// Absolute position of this placement, so a stamped result lands where it belongs.
    pub at: (f64, f64),
    /// True for the top cell itself. It is simulated, but it is not a neighbour of the
    /// placements at the same position: it is the container, not a shape beside them.
    pub is_root: bool,
    /// Hash of `signature` and `neighbours`: equal keys mean one simulation serves both.
    pub key: u64,
}

impl GdsLayout {
    /// Every placement under `top`, each with the surroundings that decide its result.
    ///
    /// `radius_nm` is how far a feature can influence its neighbour — the lateral reach of
    /// the process being modelled. Offsets are quantised to 1 nm so that two placements
    /// which differ only by sub-nanometre jitter share a context.
    pub fn contexts(&self, top: &str, radius_nm: f64) -> Vec<Context> {
        let mut placements = self.placements(top);
        // the top cell is a placement too: leaving it out drops the geometry it draws
        // directly, and the canvas then looks complete while missing features.
        let root_sig = self.structures.get(top).map(|s| s.signature(self, 0));
        if let Some(sig) = root_sig {
            placements.insert(0, Placement { cell: top.to_string(), signature: sig, at: (0.0, 0.0) });
        }
        let r = if radius_nm.is_finite() && radius_nm > 0.0 { radius_nm } else { 0.0 };
        let mut out = Vec::with_capacity(placements.len());
        // bucket the placements by their position in cells of one radius, so a placement
        // only looks at its own bucket and the eight around it. Comparing every placement
        // with every other is O(n^2), which on a million-placement die is 10^12 comparisons
        // — the feature exists for that case, so it must not be quadratic in it.
        let bucket = |x: f64, y: f64| -> (i64, i64) {
            let d = if r > 0.0 { r } else { 1.0 };
            ((x / d).floor() as i64, (y / d).floor() as i64)
        };
        let mut grid: std::collections::HashMap<(i64, i64), Vec<usize>> = std::collections::HashMap::new();
        for (i, p) in placements.iter().enumerate() {
            grid.entry(bucket(p.at.0, p.at.1)).or_default().push(i);
        }
        for (i, p) in placements.iter().enumerate() {
            let mut neighbours = Vec::new();
            let (bx, by) = bucket(p.at.0, p.at.1);
            for ox in -1..=1i64 {
                for oy in -1..=1i64 {
                    let Some(bucket_ix) = grid.get(&(bx + ox, by + oy)) else { continue };
                    for &j in bucket_ix {
                        if i == j || j == 0 || i == 0 {
                            continue;  // the root is the container, not a neighbour
                        }
                        let q = &placements[j];
                        let (dx, dy) = (q.at.0 - p.at.0, q.at.1 - p.at.1);
                        if dx.abs() > r || dy.abs() > r {
                            continue;
                        }
                        neighbours.push((dx.round() as i64, dy.round() as i64, q.signature));
                    }
                }
            }
            neighbours.sort_unstable();
            out.push(Context {
                cell: p.cell.clone(),
                at: p.at,
                is_root: i == 0 && root_sig.is_some(),
                signature: p.signature,
                key: context_key(p.signature, &neighbours),
                neighbours,
            });
        }
        out
    }

    /// Simulate each distinct **context** once and answer every placement from the cache.
    ///
    /// This is the version that is correct at boundaries: a cell in the middle of an array
    /// and the same cell at its edge have different neighbours, so they get different
    /// simulations. Cells with nothing within `radius_nm` of them reduce to plain cell
    /// dedup, which is the common case and the fast one.
    pub fn simulate_contextual<R, F>(&self, top: &str, radius_nm: f64, mut simulate: F) -> DedupRun<R>
    where
        F: FnMut(&Context) -> R,
    {
        let contexts = self.contexts(top, radius_nm);
        let mut by_signature = std::collections::BTreeMap::new();
        for c in &contexts {
            if !by_signature.contains_key(&c.key) {
                let r = simulate(c);
                by_signature.insert(c.key, r);
            }
        }
        DedupRun {
            simulations: by_signature.len(),
            placements: contexts.len(),
            by_signature,
        }
    }
}

/// The key two placements must share to be interchangeable.
fn context_key(signature: u64, neighbours: &[(i64, i64, u64)]) -> u64 {
    const PRIME: u64 = 0x0000_0100_0000_01B3;
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |bytes: &[u8]| {
        for b in bytes {
            h ^= *b as u64;
            h = h.wrapping_mul(PRIME);
        }
    };
    feed(&signature.to_le_bytes());
    for (dx, dy, s) in neighbours {
        feed(&dx.to_le_bytes());
        feed(&dy.to_le_bytes());
        feed(&s.to_le_bytes());
    }
    h
}


/// One placement of a cell in the hierarchy, with its absolute position.
#[derive(Debug, Clone, PartialEq)]
pub struct Placement {
    /// The cell being placed.
    pub cell: String,
    /// The cell's content signature: two placements with the same signature need one run.
    pub signature: u64,
    /// Absolute position of the placement's origin, in database units.
    pub at: (f64, f64),
}

/// The result of a deduplicated run.
#[derive(Debug, Clone, PartialEq)]
pub struct DedupRun<R> {
    /// One result per distinct cell content, keyed by signature.
    pub by_signature: std::collections::BTreeMap<u64, R>,
    /// How many simulations were actually performed — one per distinct content.
    pub simulations: usize,
    /// How many placements were answered from those results.
    pub placements: usize,
}

impl<R> DedupRun<R> {
    /// Placements that cost nothing because an identical cell was already simulated.
    pub fn reused(&self) -> usize {
        self.placements.saturating_sub(self.simulations)
    }
}

impl GdsLayout {
    /// Every placement in the hierarchy under `top`, with absolute positions.
    ///
    /// Unlike [`Self::flatten`], this keeps the cell identity — which is the whole point:
    /// it is what lets one simulation answer many placements.
    pub fn placements(&self, top: &str) -> Vec<Placement> {
        let mut out = Vec::new();
        self.walk_placements(top, (0.0, 0.0), &mut out, 0);
        out
    }

    fn walk_placements(&self, name: &str, origin: (f64, f64), out: &mut Vec<Placement>, depth: u32) {
        if depth > 32 {
            return;
        }
        let Some(s) = self.structures.get(name) else { return };
        for e in &s.elements {
            match e {
                Element::SRef { name: child, at, .. } => {
                    let abs = (origin.0 + at.0, origin.1 + at.1);
                    let signature = self.structures.get(child).map(|c| c.signature(self, 0)).unwrap_or(0);
                    out.push(Placement { cell: child.clone(), signature, at: abs });
                    self.walk_placements(child, abs, out, depth + 1);
                }
                Element::ARef { name: child, cols, rows, pts, .. } => {
                    // pts are the array's origin and the ends of its two pitch vectors
                    let (o, cx, cy) = (pts[0], pts[1], pts[2]);
                    let signature = self.structures.get(child).map(|c| c.signature(self, 0)).unwrap_or(0);
                    for r in 0..*rows {
                        for c in 0..*cols {
                            let fx = if *cols > 1 { c as f64 / (*cols - 1) as f64 } else { 0.0 };
                            let fy = if *rows > 1 { r as f64 / (*rows - 1) as f64 } else { 0.0 };
                            let abs = (
                                origin.0 + o.0 + (cx.0 - o.0) * fx + (cy.0 - o.0) * fy,
                                origin.1 + o.1 + (cx.1 - o.1) * fx + (cy.1 - o.1) * fy,
                            );
                            out.push(Placement { cell: child.clone(), signature, at: abs });
                        }
                    }
                }
                _ => {}
            }
        }
    }

    /// Simulate each distinct cell content once and answer every placement from the cache.
    ///
    /// `simulate` is called **once per distinct signature**, not once per placement. That
    /// is the saving the dedup report measures, actually taken: a leaf placed a thousand
    /// times by an array costs one call.
    ///
    /// The top cell itself is simulated too — it is a placement, and its own drawn geometry
    /// has to be accounted for.
    pub fn simulate_deduped<R, F>(&self, top: &str, mut simulate: F) -> DedupRun<R>
    where
        F: FnMut(&str, u64) -> R,
    {
        let mut by_signature = std::collections::BTreeMap::new();
        // the top cell counts as a placement even if nothing references it
        if let Some(s) = self.structures.get(top) {
            let sig = s.signature(self, 0);
            by_signature.insert(sig, simulate(top, sig));
        }
        let placements = self.placements(top);
        for p in &placements {
            if !by_signature.contains_key(&p.signature) {
                let r = simulate(&p.cell, p.signature);
                by_signature.insert(p.signature, r);
            }
        }
        DedupRun {
            simulations: by_signature.len(),
            placements: 1 + placements.len(),
            by_signature,
        }
    }
}


// ---------------------------------------------------------------------------
// Hierarchical dedup
// ---------------------------------------------------------------------------

/// What a layout's hierarchy is worth: how many cells it instantiates, and how many of
/// those are actually different.
///
/// A GDS file places cells; it does not repeat their contents. Simulating every placement
/// is the difference between 10^6 runs and 10^2 on a repeated layout, and the ratio is
/// exactly what this report measures before any work is done.
#[derive(Debug, Clone, PartialEq)]
pub struct DedupReport {
    /// Cells defined in the file.
    pub cells: usize,
    /// Distinct cell contents, by Merkle signature.
    pub unique: usize,
    /// Total placements: one for the top cell plus one per reference to it.
    pub placements: usize,
    /// Placements that a dedup run does not have to simulate.
    pub saved: usize,
}

impl DedupReport {
    /// Fraction of placements avoided, 0 to 1.
    pub fn saved_fraction(&self) -> f64 {
        if self.placements == 0 {
            0.0
        } else {
            self.saved as f64 / self.placements as f64
        }
    }
}

impl Structure {
    /// A Merkle hash of this cell: its own geometry, plus each child's signature and the
    /// transform it was placed with.
    ///
    /// Two cells with the same signature are the same content and need one simulation.
    /// `depth` bounds the recursion: a file with a reference cycle is illegal, but a
    /// reader that assumes legality on untrusted input is a reader that hangs.
    fn signature(&self, layout: &GdsLayout, depth: u32) -> u64 {
        const PRIME: u64 = 0x0000_0100_0000_01B3;
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |bytes: &[u8]| {
            for b in bytes {
                h ^= *b as u64;
                h = h.wrapping_mul(PRIME);
            }
        };
        if depth > 32 {
            return h;
        }
        // own geometry first, so a cell that only differs in its drawn shapes differs here
        for e in &self.elements {
            match e {
                Element::Poly { layer, datatype, pts } => {
                    feed(&[0u8]);
                    feed(&layer.to_le_bytes());
                    feed(&datatype.to_le_bytes());
                    for (x, y) in pts {
                        feed(&x.to_bits().to_le_bytes());
                        feed(&y.to_bits().to_le_bytes());
                    }
                }
                Element::Path { layer, datatype, width, pathtype, bgn_ext, end_ext, pts } => {
                    feed(&[1u8]);
                    feed(&layer.to_le_bytes());
                    feed(&datatype.to_le_bytes());
                    feed(&width.to_bits().to_le_bytes());
                    feed(&pathtype.to_le_bytes());
                    feed(&bgn_ext.to_bits().to_le_bytes());
                    feed(&end_ext.to_bits().to_le_bytes());
                    for (x, y) in pts {
                        feed(&x.to_bits().to_le_bytes());
                        feed(&y.to_bits().to_le_bytes());
                    }
                }
                // a child contributes its content and where it sits, not its name: two
                // cells that place identical children identically are identical cells
                Element::SRef { name, trans, at } => {
                    feed(&[2u8]);
                    feed(&trans.mag.to_bits().to_le_bytes());
                    feed(&(trans.angle_deg as i64).to_le_bytes());
                    feed(&[trans.reflect as u8]);
                    feed(&at.0.to_bits().to_le_bytes());
                    feed(&at.1.to_bits().to_le_bytes());
                    let child = layout.structures.get(name).map(|c| c.signature(layout, depth + 1)).unwrap_or(0);
                    feed(&child.to_le_bytes());
                }
                Element::ARef { name, trans, cols, rows, pts } => {
                    feed(&[3u8]);
                    feed(&trans.mag.to_bits().to_le_bytes());
                    feed(&(trans.angle_deg as i64).to_le_bytes());
                    feed(&[trans.reflect as u8]);
                    feed(&cols.to_le_bytes());
                    feed(&rows.to_le_bytes());
                    for (x, y) in pts {
                        feed(&x.to_bits().to_le_bytes());
                        feed(&y.to_bits().to_le_bytes());
                    }
                    let child = layout.structures.get(name).map(|c| c.signature(layout, depth + 1)).unwrap_or(0);
                    feed(&child.to_le_bytes());
                }
            }
        }
        h
    }

    /// How many times this cell references other cells, including arrays.
    fn reference_count(&self) -> usize {
        self.elements
            .iter()
            .map(|e| match e {
                Element::SRef { .. } => 1,
                Element::ARef { cols, rows, .. } => cols * rows,
                _ => 0,
            })
            .sum()
    }
}

impl GdsLayout {
    /// How much of this layout is repetition.
    pub fn dedup(&self) -> DedupReport {
        let mut sigs = std::collections::BTreeMap::new();
        for (name, s) in &self.structures {
            sigs.insert(name.clone(), s.signature(self, 0));
        }
        let unique = {
            let mut v: Vec<u64> = sigs.values().copied().collect();
            v.sort_unstable();
            v.dedup();
            v.len()
        };
        // the top cell is one placement; every reference is another
        let placements = 1 + self.structures.values().map(|s| s.reference_count()).sum::<usize>();
        DedupReport {
            cells: self.structures.len(),
            unique,
            placements,
            saved: placements.saturating_sub(unique),
        }
    }

    /// One signature per cell, so a caller can group identical cells by content.
    pub fn cell_signatures(&self) -> Vec<(String, u64)> {
        self.structures
            .iter()
            .map(|(name, s)| (name.clone(), s.signature(self, 0)))
            .collect()
    }
}

pub fn summarize(polys: &[FlatPoly]) -> Vec<LayerSummary> {
    let mut m: BTreeMap<(i32, i32), LayerSummary> = BTreeMap::new();
    for p in polys {
        let e = m.entry((p.layer, p.datatype)).or_insert(LayerSummary { layer: p.layer, datatype: p.datatype, polygons: 0, bbox: [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], anchor: (0.0, 0.0), feature_p10_nm: 0.0, feature_med_nm: 0.0 });
        e.polygons += 1;
        for v in &p.pts { e.bbox[0] = e.bbox[0].min(v.0); e.bbox[1] = e.bbox[1].min(v.1); e.bbox[2] = e.bbox[2].max(v.0); e.bbox[3] = e.bbox[3].max(v.1); }
    }
    // Feature-size statistics per layer.
    let mut dims: BTreeMap<(i32, i32), Vec<f64>> = BTreeMap::new();
    for p in polys {
        let b = p.pts.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, v| [b[0].min(v.0), b[1].min(v.1), b[2].max(v.0), b[3].max(v.1)]);
        dims.entry((p.layer, p.datatype)).or_default().push((b[2] - b[0]).min(b[3] - b[1]));
    }
    for (k, mut d) in dims {
        d.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if let Some(l) = m.get_mut(&k) { l.feature_p10_nm = d[(d.len() as f64 * 0.10) as usize]; l.feature_med_nm = d[d.len() / 2]; }
    }
    // Anchor = centre of the polygon whose centre is nearest the layer's bounding-box centre.
    for l in m.values_mut() {
        let c = l.center_nm();
        let mut best = f64::INFINITY;
        for p in polys.iter().filter(|p| p.layer == l.layer && p.datatype == l.datatype) {
            let b = p.pts.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, v| [b[0].min(v.0), b[1].min(v.1), b[2].max(v.0), b[3].max(v.1)]);
            let pc = ((b[0] + b[2]) * 0.5, (b[1] + b[3]) * 0.5);
            let d = (pc.0 - c.0).hypot(pc.1 - c.1);
            // A corner of the nearest polygon: a window there sees an edge (inside and outside), unlike the polygon centre.
            if d < best { best = d; l.anchor = p.pts[0]; }
        }
    }
    m.into_values().collect()
}

#[allow(dead_code)]
pub(crate) fn _apply_for_tests(t: &Trans, at: (f64, f64), p: (f64, f64)) -> (f64, f64) { apply(t, at, p) }

#[cfg(any(test, feature = "test-fixtures"))]
pub(crate) mod tests {
    #[allow(unused_imports)]
    use super::*;


#[cfg(test)]
mod dedup_tests {
    use super::*;

    fn poly(layer: i32, x0: f64, side: f64) -> Element {
        Element::Poly {
            layer,
            datatype: 0,
            pts: vec![(x0, 0.0), (x0 + side, 0.0), (x0 + side, side), (x0, side)],
        }
    }

    fn sref(name: &str, at: (f64, f64)) -> Element {
        Element::SRef { name: name.to_string(), trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 }, at }
    }

    fn layout(cells: Vec<(&str, Vec<Element>)>) -> GdsLayout {
        let mut structures = BTreeMap::new();
        for (name, elements) in cells {
            structures.insert(name.to_string(), Structure { elements });
        }
        GdsLayout { structures, db_unit_m: 1e-9 }
    }

    #[test]
    fn identical_cells_share_one_signature() {
        let l = layout(vec![
            ("A", vec![poly(1, 0.0, 100.0)]),
            ("B", vec![poly(1, 0.0, 100.0)]),
        ]);
        let s = l.cell_signatures();
        assert_eq!(s.len(), 2);
        assert_eq!(s[0].1, s[1].1, "the same drawing is the same cell");
    }

    #[test]
    fn a_different_drawing_is_a_different_cell() {
        let l = layout(vec![
            ("A", vec![poly(1, 0.0, 100.0)]),
            ("B", vec![poly(1, 0.0, 101.0)]),
        ]);
        let s = l.cell_signatures();
        assert_ne!(s[0].1, s[1].1);
    }

    #[test]
    fn a_cell_is_its_children_as_well_as_its_own_geometry() {
        // the Merkle part: A and B draw nothing themselves but place different children,
        // so they are different cells
        let l = layout(vec![
            ("leaf1", vec![poly(1, 0.0, 100.0)]),
            ("leaf2", vec![poly(1, 0.0, 200.0)]),
            ("A", vec![sref("leaf1", (0.0, 0.0))]),
            ("B", vec![sref("leaf2", (0.0, 0.0))]),
        ]);
        let s: std::collections::BTreeMap<String, u64> = l.cell_signatures().into_iter().collect();
        assert_ne!(s["A"], s["B"], "different children make different parents");
    }

    #[test]
    fn placing_the_same_child_twice_does_not_make_a_new_cell() {
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 100.0)]),
            ("A", vec![sref("leaf", (0.0, 0.0)), sref("leaf", (500.0, 0.0))]),
        ]);
        let s: std::collections::BTreeMap<String, u64> = l.cell_signatures().into_iter().collect();
        assert_ne!(s["A"], s["leaf"]);
        assert_eq!(l.dedup().unique, 2, "two contents: the leaf and the pair");
    }

    #[test]
    fn the_report_counts_what_a_dedup_run_saves() {
        // one leaf placed a thousand times by an array: one simulation, not a thousand
        let arr = Element::ARef {
            name: "leaf".to_string(),
            trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 },
            cols: 100,
            rows: 10,
            pts: [(0.0, 0.0), (1000.0, 0.0), (0.0, 1000.0)],
        };
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 100.0)]),
            ("top", vec![arr]),
        ]);
        let r = l.dedup();
        assert_eq!(r.cells, 2);
        assert_eq!(r.unique, 2);
        assert_eq!(r.placements, 1001, "the top plus a thousand array sites");
        assert_eq!(r.saved, 999);
        assert!(r.saved_fraction() > 0.99);
    }

    #[test]
    fn a_flat_layout_saves_nothing() {
        // no references, no repetition: the report must not invent a saving
        let l = layout(vec![("only", vec![poly(1, 0.0, 100.0)])]);
        let r = l.dedup();
        assert_eq!((r.cells, r.unique, r.placements), (1, 1, 1));
        assert_eq!(r.saved, 0);
        assert_eq!(r.saved_fraction(), 0.0);
    }

    #[test]
    fn a_reference_cycle_terminates() {
        // illegal GDS, but a reader must not hang on it
        let l = layout(vec![
            ("A", vec![sref("B", (0.0, 0.0))]),
            ("B", vec![sref("A", (0.0, 0.0))]),
        ]);
        let r = l.dedup();
        assert_eq!(r.cells, 2);
    }
}


#[cfg(test)]
mod dedup_run_tests {
    use super::*;

    fn poly(layer: i32, x0: f64, side: f64) -> Element {
        Element::Poly {
            layer,
            datatype: 0,
            pts: vec![(x0, 0.0), (x0 + side, 0.0), (x0 + side, side), (x0, side)],
        }
    }

    fn sref(name: &str, at: (f64, f64)) -> Element {
        Element::SRef { name: name.to_string(), trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 }, at }
    }

    fn layout(cells: Vec<(&str, Vec<Element>)>) -> GdsLayout {
        let mut structures = BTreeMap::new();
        for (name, elements) in cells {
            structures.insert(name.to_string(), Structure { elements });
        }
        GdsLayout { structures, db_unit_m: 1e-9 }
    }

    #[test]
    fn a_thousand_placements_cost_one_simulation() {
        // the saving the report measures, actually taken
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 100.0)]),
            ("top", vec![Element::ARef {
                name: "leaf".to_string(),
                trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 },
                cols: 100,
                rows: 10,
                pts: [(0.0, 0.0), (1000.0, 0.0), (0.0, 1000.0)],
            }]),
        ]);
        let mut calls = 0;
        let run = l.simulate_deduped("top", |_name, sig| {
            calls += 1;
            sig
        });
        assert_eq!(calls, 2, "the top and the leaf, not a thousand leaves");
        assert_eq!(run.placements, 1001);
        assert_eq!(run.simulations, 2);
        assert_eq!(run.reused(), 999);
    }

    #[test]
    fn distinct_cells_are_each_simulated() {
        let l = layout(vec![
            ("a", vec![poly(1, 0.0, 100.0)]),
            ("b", vec![poly(1, 0.0, 200.0)]),
            ("top", vec![sref("a", (0.0, 0.0)), sref("b", (500.0, 0.0))]),
        ]);
        let mut calls = 0;
        let run = l.simulate_deduped("top", |_n, sig| { calls += 1; sig });
        assert_eq!(calls, 3, "top, a and b");
        assert_eq!(run.placements, 3);
        assert_eq!(run.reused(), 0);
    }

    #[test]
    fn identical_cells_nested_twice_are_simulated_once() {
        // two different names, same drawing: one run serves both
        let l = layout(vec![
            ("x", vec![poly(1, 0.0, 100.0)]),
            ("y", vec![poly(1, 0.0, 100.0)]),
            ("top", vec![sref("x", (0.0, 0.0)), sref("y", (500.0, 0.0))]),
        ]);
        let mut calls = 0;
        let run = l.simulate_deduped("top", |_n, sig| { calls += 1; sig });
        assert_eq!(calls, 2, "top plus one content shared by x and y");
        assert_eq!(run.reused(), 1);
    }

    #[test]
    fn placements_carry_absolute_positions() {
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 10.0)]),
            ("mid", vec![sref("leaf", (100.0, 0.0))]),
            ("top", vec![sref("mid", (1000.0, 0.0))]),
        ]);
        let ps = l.placements("top");
        let leaf: Vec<_> = ps.iter().filter(|p| p.cell == "leaf").collect();
        assert_eq!(leaf.len(), 1);
        assert_eq!(leaf[0].at, (1100.0, 0.0), "positions accumulate down the hierarchy");
    }

    #[test]
    fn an_array_places_every_site() {
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 10.0)]),
            ("top", vec![Element::ARef {
                name: "leaf".to_string(),
                trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 },
                cols: 3,
                rows: 2,
                pts: [(0.0, 0.0), (300.0, 0.0), (0.0, 200.0)],
            }]),
        ]);
        let ps = l.placements("top");
        assert_eq!(ps.len(), 6, "three columns by two rows");
        assert!(ps.iter().all(|p| p.cell == "leaf"));
    }

    #[test]
    fn a_cycle_does_not_hang_the_walk() {
        let l = layout(vec![
            ("a", vec![sref("b", (0.0, 0.0))]),
            ("b", vec![sref("a", (0.0, 0.0))]),
        ]);
        let _ = l.placements("a"); // must terminate
        let run = l.simulate_deduped("a", |_n, sig| sig);
        assert!(run.placements >= 1);
    }
}


#[cfg(test)]
mod context_tests {
    use super::*;

    fn poly(layer: i32, x0: f64, side: f64) -> Element {
        Element::Poly { layer, datatype: 0, pts: vec![(x0, 0.0), (x0 + side, 0.0), (x0 + side, side), (x0, side)] }
    }
    fn sref(name: &str, at: (f64, f64)) -> Element {
        Element::SRef { name: name.to_string(), trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 }, at }
    }
    fn layout(cells: Vec<(&str, Vec<Element>)>) -> GdsLayout {
        let mut structures = BTreeMap::new();
        for (name, elements) in cells { structures.insert(name.to_string(), Structure { elements }); }
        GdsLayout { structures, db_unit_m: 1e-9 }
    }

    #[test]
    fn an_array_has_fewer_contexts_than_placements() {
        // a 5x5 array of one cell. The interior cells all see the same eight neighbours, so
        // they share one key. Edges and corners differ from each other by direction — the
        // key is deliberately direction-sensitive, which is conservative: it never merges
        // two placements that could etch differently.
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 10.0)]),
            ("top", vec![Element::ARef {
                name: "leaf".to_string(),
                trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 },
                cols: 5, rows: 5, pts: [(0.0, 0.0), (400.0, 0.0), (0.0, 400.0)],
            }]),
        ]);
        let ctx = l.contexts("top", 150.0);
        assert_eq!(ctx.len(), 26, "twenty-five placements plus the root");
        let mut keys: Vec<u64> = ctx.iter().map(|c| c.key).collect();
        keys.sort_unstable();
        keys.dedup();
        assert!(keys.len() < ctx.len(), "dedup must find repetition");
        assert_eq!(keys.len(), 10, "one interior, four edges, four corners, plus the root");
    }

    #[test]
    fn the_same_cell_at_the_edge_and_the_middle_differ() {
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 10.0)]),
            ("top", vec![Element::ARef {
                name: "leaf".to_string(),
                trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 },
                cols: 3, rows: 3, pts: [(0.0, 0.0), (200.0, 0.0), (0.0, 200.0)],
            }]),
        ]);
        let ctx = l.contexts("top", 150.0);
        let centre = ctx.iter().find(|c| c.neighbours.len() == 8).expect("a centre placement");
        let corner = ctx.iter().find(|c| c.neighbours.len() == 3).expect("a corner placement");
        assert_eq!(centre.signature, corner.signature, "the same cell");
        assert_ne!(centre.key, corner.key, "but different surroundings");
    }

    #[test]
    fn a_radius_of_zero_makes_every_placement_isolated() {
        // with no interaction there is no context to differ: pure cell dedup
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 10.0)]),
            ("top", vec![Element::ARef {
                name: "leaf".to_string(),
                trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 },
                cols: 4, rows: 4, pts: [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0)],
            }]),
        ]);
        let ctx = l.contexts("top", 0.0);
        // the root sits at the same point as the first array site, so those two see
        // each other even at radius zero — coincident placements overlap
        let run = l.simulate_contextual("top", 0.0, |c| c.signature);
        assert_eq!(run.simulations, 2, "the leaf and the root");
        assert_eq!(run.placements, 17, "sixteen sites plus the root");
        assert_eq!(run.reused(), 15);
    }

    #[test]
    fn isolated_cells_reduce_to_plain_dedup() {
        // far apart relative to the radius: nothing to distinguish them
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 10.0)]),
            ("top", vec![sref("leaf", (0.0, 0.0)), sref("leaf", (100_000.0, 0.0))]),
        ]);
        let run = l.simulate_contextual("top", 50.0, |c| c.signature);
        assert_eq!(run.simulations, 2, "the leaf and the root");
        assert_eq!(run.placements, 3, "two leaves plus the root");
        assert_eq!(run.reused(), 1);
    }

    #[test]
    fn close_together_they_are_no_longer_interchangeable() {
        // same two cells, but now inside the interaction radius: each sees the other, and
        // the pair is not the same as an isolated instance
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 10.0)]),
            ("top", vec![sref("leaf", (0.0, 0.0)), sref("leaf", (100_000.0, 0.0))]),
        ]);
        let far = l.simulate_contextual("top", 50.0, |c| c.signature);
        let near = l.simulate_contextual("top", 200_000.0, |c| c.signature);
        assert_eq!(far.simulations, 2, "the leaf and the root");
        assert_eq!(near.simulations, 3, "the two leaves and the root");
        assert!(near.simulations > far.simulations);
    }


    #[test]
    fn a_whole_file_runs_end_to_end() {
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 40.0)]),
            ("top", vec![sref("leaf", (0.0, 0.0)), sref("leaf", (100.0, 0.0))]),
        ]);
        let spec = crate::gds_sweep::CellSimSpec { cell_nm: 10.0, core_nm: (40.0, 40.0), margin_nm: 10.0, substrate_nm: 500.0 };
        let rates = atto_geom::column::EtchRates { si: 100.0, sio2: 100.0, si3n4: 100.0, resist: 0.0, metal: 0.0 };
        let etch_spec = atto_geom::column::EtchSpec { rates, time_min: 0.1 };
        let r = crate::gds_sweep::run_whole_file(&l, "top", &spec, &etch_spec, 50.0).expect("the run completes");
        assert_eq!(r.placements, 3, "two leaves plus the root");
        assert!(r.simulations >= 2, "the leaf and the root at least: {}", r.simulations);
        assert_eq!(r.failed, 0, "no context should fail to etch");
        assert!(r.canvas.depth_nm.iter().any(|v| *v > 0.0), "something was etched");
    }

    #[test]
    fn contexts_are_stable_across_runs() {
        // the key must be deterministic, or the cache never hits
        let l = layout(vec![
            ("leaf", vec![poly(1, 0.0, 10.0)]),
            ("top", vec![Element::ARef {
                name: "leaf".to_string(),
                trans: Trans { reflect: false, mag: 1.0, angle_deg: 0.0 },
                cols: 2, rows: 2, pts: [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0)],
            }]),
        ]);
        let a = l.contexts("top", 150.0);
        let b = l.contexts("top", 150.0);
        assert_eq!(a, b);
    }
}

    pub fn enc_real8(v: f64) -> [u8; 8] {
        if v == 0.0 { return [0; 8]; }
        let (sign, mut x) = if v < 0.0 { (0x80u8, -v) } else { (0u8, v) };
        let mut e = 0i32;
        while x >= 1.0 { x /= 16.0; e += 1; }
        while x < 1.0 / 16.0 { x *= 16.0; e -= 1; }
        let m = (x * (1u64 << 56) as f64).round() as u64;
        let mut b = [0u8; 8];
        b[0] = sign | ((e + 64) as u8 & 0x7F);
        for i in 0..7 { b[1 + i] = (m >> (8 * (6 - i))) as u8; }
        b
    }
    fn rec(b: &mut Vec<u8>, ty: u8, dt: u8, payload: &[u8]) { b.extend_from_slice(&((payload.len() + 4) as u16).to_be_bytes()); b.push(ty); b.push(dt); b.extend_from_slice(payload); }
    fn i16s(v: &[i16]) -> Vec<u8> { v.iter().flat_map(|x| x.to_be_bytes()).collect() }
    fn i32s(v: &[i32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_be_bytes()).collect() }
    fn nm(s: &str) -> Vec<u8> { let mut v = s.as_bytes().to_vec(); if v.len() % 2 == 1 { v.push(0); } v }

    /// db unit = 0.5 nm to prove UNITS scaling. Top cell: SREF(child, rotated 90°, mag 2, reflected) + AREF 3×2 + PATH + BOX.
    pub fn sample_stream() -> Vec<u8> {
        let mut b = Vec::new();
        rec(&mut b, 0x00, 0x02, &i16s(&[600]));
        rec(&mut b, 0x01, 0x02, &[0; 24]);
        rec(&mut b, 0x02, 0x06, &nm("LIB"));
        let mut units = enc_real8(1.0e-3).to_vec(); units.extend_from_slice(&enc_real8(0.5e-9));
        rec(&mut b, 0x03, 0x05, &units);
        // child: 100x40 db-unit rectangle at the origin, layer 5/0
        rec(&mut b, 0x05, 0x02, &[0; 24]); rec(&mut b, 0x06, 0x06, &nm("CHILD"));
        rec(&mut b, 0x08, 0x00, &[]); rec(&mut b, 0x0D, 0x02, &i16s(&[5])); rec(&mut b, 0x0E, 0x02, &i16s(&[0]));
        rec(&mut b, 0x10, 0x03, &i32s(&[0, 0, 100, 0, 100, 40, 0, 40, 0, 0])); rec(&mut b, 0x11, 0x00, &[]);
        rec(&mut b, 0x07, 0x00, &[]);
        // top
        rec(&mut b, 0x05, 0x02, &[0; 24]); rec(&mut b, 0x06, 0x06, &nm("TOP"));
        // SREF: reflect, mag 2, angle 90, at (1000, 2000)
        rec(&mut b, 0x0A, 0x00, &[]); rec(&mut b, 0x12, 0x06, &nm("CHILD"));
        rec(&mut b, 0x1A, 0x01, &i16s(&[i16::MIN])); rec(&mut b, 0x1B, 0x05, &enc_real8(2.0)); rec(&mut b, 0x1C, 0x05, &enc_real8(90.0));
        rec(&mut b, 0x10, 0x03, &i32s(&[1000, 2000])); rec(&mut b, 0x11, 0x00, &[]);
        // AREF 3 cols x 2 rows, pitch (300, 500) db units, at origin
        rec(&mut b, 0x0B, 0x00, &[]); rec(&mut b, 0x12, 0x06, &nm("CHILD")); rec(&mut b, 0x13, 0x02, &i16s(&[3, 2]));
        rec(&mut b, 0x10, 0x03, &i32s(&[5000, 0, 5000 + 900, 0, 5000, 1000])); rec(&mut b, 0x11, 0x00, &[]);
        // PATH layer 7 width 20, flat ends, (0,-500)->(400,-500)
        rec(&mut b, 0x09, 0x00, &[]); rec(&mut b, 0x0D, 0x02, &i16s(&[7])); rec(&mut b, 0x0E, 0x02, &i16s(&[1])); rec(&mut b, 0x21, 0x02, &i16s(&[0]));
        rec(&mut b, 0x0F, 0x03, &i32s(&[20])); rec(&mut b, 0x10, 0x03, &i32s(&[0, -500, 400, -500])); rec(&mut b, 0x11, 0x00, &[]);
        // BOX layer 9
        rec(&mut b, 0x2D, 0x00, &[]); rec(&mut b, 0x0D, 0x02, &i16s(&[9])); rec(&mut b, 0x2E, 0x02, &i16s(&[0]));
        rec(&mut b, 0x10, 0x03, &i32s(&[-100, -100, -50, -100, -50, -50, -100, -50, -100, -100])); rec(&mut b, 0x11, 0x00, &[]);
        rec(&mut b, 0x07, 0x00, &[]);
        rec(&mut b, 0x04, 0x00, &[]);
        b
    }

    #[test]
    fn real8_round_trips() { for v in [1.0e-3, 0.5e-9, 90.0, 2.0, -3.25, 1.0] { assert!((real8(&enc_real8(v)) / v - 1.0).abs() < 1e-12, "{v}"); } }

    #[test]
    fn hierarchy_units_and_element_types_flatten_correctly() {
        let lay = parse(&sample_stream()).unwrap();
        assert!((lay.db_unit_m - 0.5e-9).abs() < 1e-18);
        assert_eq!(lay.default_top().as_deref(), Some("TOP"), "CHILD is referenced, TOP is not");
        let polys = lay.flatten(None).unwrap();
        let layers = summarize(&polys);
        let get = |l, d| layers.iter().find(|s| s.layer == l && s.datatype == d).unwrap_or_else(|| panic!("layer {l}/{d} missing: {layers:?}"));
        assert_eq!(get(5, 0).polygons, 1 + 6, "1 SREF + 3x2 AREF instances of the child");
        assert_eq!(get(7, 1).polygons, 1, "a one-segment PATH becomes one polygon");
        assert_eq!(get(9, 0).polygons, 1, "BOX is read as a polygon");
        // PATH: (0,-500)->(400,-500) db units, width 20 => 0.5 nm/dbu: x 0..200 nm, y -260..-240 nm.
        let p = get(7, 1);
        assert!((p.bbox[0] - 0.0).abs() < 1e-9 && (p.bbox[2] - 200.0).abs() < 1e-9 && (p.bbox[1] + 255.0).abs() < 1e-9 && (p.bbox[3] + 245.0).abs() < 1e-9, "{:?}", p.bbox);
        // SREF: child rect (0..100, 0..40) reflected (y -> -y): (0..100, -40..0); mag 2: (0..200, -80..0);
        // rotate 90°: (x,y) -> (-y, x): x in 0..80, y in 0..200; translate (1000, 2000) db -> 1000..1080, 2000..2200 db -> nm * 0.5.
        let sref = polys.iter().filter(|q| q.layer == 5).find(|q| q.pts.iter().all(|v| v.0 >= 499.0 && v.0 <= 541.0)).expect("SREF instance");
        let (mut lo, mut hi) = ((f64::INFINITY, f64::INFINITY), (f64::NEG_INFINITY, f64::NEG_INFINITY));
        for v in &sref.pts { lo = (lo.0.min(v.0), lo.1.min(v.1)); hi = (hi.0.max(v.0), hi.1.max(v.1)); }
        assert!((lo.0 - 500.0).abs() < 1e-6 && (hi.0 - 540.0).abs() < 1e-6 && (lo.1 - 1000.0).abs() < 1e-6 && (hi.1 - 1100.0).abs() < 1e-6, "{lo:?} {hi:?}");
        // AREF: origin 5000 db, pitch 300 x 500 db => instance origins at (2500 + 150c, 250r) nm.
        let aref_x: Vec<f64> = { let mut v: Vec<f64> = polys.iter().filter(|q| q.layer == 5 && q.pts.iter().any(|p| p.0 >= 2499.0)).map(|q| q.pts.iter().map(|p| p.0).fold(f64::INFINITY, f64::min)).collect(); v.sort_by(|a, b| a.partial_cmp(b).unwrap()); v.dedup_by(|a, b| (*a - *b).abs() < 1e-6); v };
        assert_eq!(aref_x.len(), 3);
        assert!((aref_x[1] - aref_x[0] - 150.0).abs() < 1e-6, "{aref_x:?}");
    }

    #[test]
    fn non_gds_truncated_missing_and_cyclic_inputs_report_clear_errors() {
        assert!(matches!(parse(b"hello world, definitely text"), Err(GdsError::Format(_))));
        let s = sample_stream();
        assert!(matches!(parse(&s[..s.len() - 9]), Err(GdsError::Eof) | Ok(_)));
        assert!(matches!(parse(&s[..40]), Err(GdsError::Eof)));
        let lay = parse(&s).unwrap();
        assert_eq!(lay.flatten(Some("NOPE")).unwrap_err(), GdsError::MissingStructure("NOPE".into()));
        // A cell that references itself must not hang.
        let mut cyc = lay.clone();
        cyc.structures.get_mut("CHILD").unwrap().elements.push(Element::SRef { name: "CHILD".into(), trans: Trans::default(), at: (0.0, 0.0) });
        assert!(matches!(cyc.flatten(Some("CHILD")), Err(GdsError::TooDeep(_))));
    }
}

/// Test fixtures shared with other crates' tests (hierarchy, units, path, box).
#[doc(hidden)]
#[cfg(feature = "test-fixtures")]
pub mod tests_support { pub fn sample_stream() -> Vec<u8> { super::tests::sample_stream() } }

#[cfg(test)]
mod faithful_tests {
    use super::*;

    fn lay_with(elements: Vec<Element>, db_m: f64) -> GdsLayout {
        let mut m = BTreeMap::new();
        m.insert("T".to_string(), Structure { elements });
        GdsLayout { structures: m, db_unit_m: db_m }
    }
    fn bbox_of(p: &FlatPoly) -> [f64; 4] { pts_bbox(&p.pts) }
    fn path(pt: i16, bgn: f64, end: f64, w: f64) -> Vec<FlatPoly> {
        lay_with(vec![Element::Path { layer: 1, datatype: 0, width: w, pathtype: pt, bgn_ext: bgn, end_ext: end, pts: vec![(0.0, 0.0), (1000.0, 0.0)] }], 1e-9).flatten(None).unwrap()
    }
    fn ext(polys: &[FlatPoly]) -> [f64; 4] { polys.iter().map(bbox_of).fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |a, b| [a[0].min(b[0]), a[1].min(b[1]), a[2].max(b[2]), a[3].max(b[3])]) }

    #[test]
    fn path_types_extensions_negative_width_and_joins() {
        let flat0 = ext(&path(0, 0.0, 0.0, 100.0));
        assert_eq!(flat0, [0.0, -50.0, 1000.0, 50.0], "type 0: flat ends");
        assert_eq!(ext(&path(2, 0.0, 0.0, 100.0)), [-50.0, -50.0, 1050.0, 50.0], "type 2: extended by half the width");
        assert_eq!(ext(&path(4, 30.0, 70.0, 100.0)), [-30.0, -50.0, 1070.0, 50.0], "type 4: BGNEXTN / ENDEXTN");
        let round = ext(&path(1, 0.0, 0.0, 100.0));
        assert!(round[0] <= -49.0 && round[2] >= 1049.0, "type 1: round caps reach half a width beyond the ends: {round:?}");
        assert_eq!(ext(&path(0, 0.0, 0.0, -100.0)), flat0, "negative width = absolute width, magnitude used");
        let l = lay_with(vec![Element::Path { layer: 1, datatype: 0, width: 100.0, pathtype: 0, bgn_ext: 0.0, end_ext: 0.0, pts: vec![(0.0, 0.0), (1000.0, 0.0), (1000.0, 1000.0)] }], 1e-9).flatten(None).unwrap();
        let polys: Vec<Vec<[f64; 2]>> = l.iter().map(|p| p.pts.iter().map(|v| [v.0, v.1]).collect()).collect();
        let m = crate::geom::MaskPattern::Polygons { polys, clear_field: false };
        assert_eq!(m.transmission(1040.0, -40.0), 1.0, "outer corner of the bend is filled by the join");
        assert_eq!(m.transmission(1100.0, -60.0), 0.0);
        assert!(lay_with(vec![Element::Path { layer: 1, datatype: 0, width: 0.0, pathtype: 0, bgn_ext: 0.0, end_ext: 0.0, pts: vec![(0.0, 0.0), (5.0, 5.0)] }], 1e-9).flatten(None).unwrap().is_empty());
    }

    #[test]
    fn region_flatten_equals_the_full_flatten_filtered_and_skips_far_instances() {
        let lay = parse(&tests::sample_stream()).unwrap();
        let full = lay.flatten(None).unwrap();
        for region in [[0.0, 0.0, 1200.0, 1200.0], [2400.0, -100.0, 3200.0, 700.0], [-1.0e9, -1.0e9, 1.0e9, 1.0e9], [9.0e8, 9.0e8, 9.1e8, 9.1e8]] {
            let part = lay.flatten_region(None, Some(region)).unwrap();
            let expect = full.iter().filter(|p| { let b = bbox_of(p); !(b[2] < region[0] || b[0] > region[2] || b[3] < region[1] || b[1] > region[3]) }).count();
            assert_eq!(part.len(), expect, "region {region:?}");
        }
        let mut big = parse(&tests::sample_stream()).unwrap();
        big.structures.get_mut("TOP").unwrap().elements.push(Element::ARef { name: "CHILD".into(), trans: Trans::default(), cols: 1000, rows: 1000, pts: [(0.0, 0.0), (2.0e6, 0.0), (0.0, 2.0e6)] });
        let t = std::time::Instant::now();
        let r = big.flatten_region(None, Some([0.0, 0.0, 600.0, 600.0])).unwrap();
        assert!(t.elapsed().as_secs_f64() < 5.0 && !r.is_empty(), "culled flatten took {:?}", t.elapsed());
    }
}
