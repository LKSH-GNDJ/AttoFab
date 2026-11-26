//! Translates the SVO wafer state into flat buffers a wgpu compute/fragment
//! raymarching shader can consume. This module is deliberately GPU-API
//! agnostic (no `wgpu` dependency in `sim_core`) — it produces plain
//! `Vec<u8>` / `Vec<f32>` buffers; the platform shells (`app_desktop`,
//! `app_web`) own the actual `wgpu::Device`/`Buffer` creation and upload.

use crate::wafer::{Bounds, Material, SvoNode, WaferState};
use bytemuck::{Pod, Zeroable};
use std::collections::VecDeque;

/// One flattened node record, matching the layout consumed by
/// `assets/shaders/wafer_render.wgsl`'s storage buffer.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct GpuSvoNode {
    /// Bounds min (x, y, z) in nm; `_pad0` supplies WGSL's 16-byte alignment.
    pub min: [f32; 3],
    pub _pad0: f32,
    /// Bounds max (x, y, z) in nm. WGSL places this at byte offset 16.
    pub max: [f32; 3],
    /// Material id (see `material_to_id`); -1 if this is an internal node.
    pub material_id: i32,
    /// Index of first child in the flattened buffer, or -1 for a leaf.
    pub first_child: i32,
    /// Explicit padding to match the 48-byte WGSL storage-buffer stride.
    pub _padding: [i32; 3],
}

/// Flatten the SVO into a linear buffer suitable for GPU upload, via a
/// breadth-first walk so parent indices are always emitted before their
/// children's own subtrees are fully resolved (each node still stores a
/// direct `first_child` pointer once its children are appended).
pub fn flatten_for_gpu(wafer: &WaferState) -> Vec<GpuSvoNode> {
    flatten_for_gpu_relative(wafer, (0.0, 0.0, 0.0))
}

/// Like [`flatten_for_gpu`] but with every node bound expressed relative to
/// `origin` (wafer-local nm) in f64 *before* the f32 cast. Rendering a tile
/// 100 mm from the wafer centre would otherwise lose ~8 nm of precision.
pub fn flatten_for_gpu_relative(wafer: &WaferState, origin: (f64, f64, f64)) -> Vec<GpuSvoNode> {
    let rel = |b: Bounds| Bounds {
        min: (b.min.0 - origin.0, b.min.1 - origin.1, b.min.2 - origin.2),
        max: (b.max.0 - origin.0, b.max.1 - origin.1, b.max.2 - origin.2),
    };
    let root_bounds = rel(wafer.bounds);
    let mut out = vec![gpu_record(&wafer.root, root_bounds)];
    let mut queue: VecDeque<(usize, &SvoNode, Bounds)> = VecDeque::new();
    queue.push_back((0, &wafer.root, root_bounds));

    // Breadth-first flattening guarantees that every node's eight direct
    // children are contiguous at `first_child .. first_child + 8`, exactly
    // as the WGSL traversal expects. A recursive DFS layout does not satisfy
    // that contract because grandchildren become interleaved with siblings.
    while let Some((index, node, bounds)) = queue.pop_front() {
        if let SvoNode::Internal(children) = node {
            let first_child = out.len() as i32;
            out[index].first_child = first_child;
            for (child_index, child) in children.iter().enumerate() {
                let child_bounds = octant_bounds(bounds, child_index);
                let flat_index = out.len();
                out.push(gpu_record(child, child_bounds));
                queue.push_back((flat_index, child, child_bounds));
            }
        }
    }
    out
}

fn gpu_record(node: &SvoNode, bounds: Bounds) -> GpuSvoNode {
    let material_id = match node {
        SvoNode::Leaf(material) => material_to_id(*material),
        SvoNode::Internal(_) => -1,
    };
    GpuSvoNode {
        min: [bounds.min.0 as f32, bounds.min.1 as f32, bounds.min.2 as f32],
        _pad0: 0.0,
        max: [bounds.max.0 as f32, bounds.max.1 as f32, bounds.max.2 as f32],
        material_id,
        first_child: -1,
        _padding: [0; 3],
    }
}

fn octant_bounds(b: Bounds, i: usize) -> Bounds {
    let c = b.center();
    let (x_lo, x_hi) = if i & 1 == 0 { (b.min.0, c.0) } else { (c.0, b.max.0) };
    let (y_lo, y_hi) = if i & 2 == 0 { (b.min.1, c.1) } else { (c.1, b.max.1) };
    let (z_lo, z_hi) = if i & 4 == 0 { (b.min.2, c.2) } else { (c.2, b.max.2) };
    Bounds { min: (x_lo, y_lo, z_lo), max: (x_hi, y_hi, z_hi) }
}

/// Number of palette slots shared by the shader and the egui legend.
pub const PALETTE_SLOTS: usize = 16;
/// First palette slot used by dopant concentration buckets.
pub const DOPANT_ID_BASE: i32 = 7;

/// Maps a Material to the palette slot the WGSL shader indexes. Dopant
/// concentration buckets 0..=8 occupy slots 7..=15.
pub fn material_to_id(material: Material) -> i32 {
    match material {
        Material::Air => 0,
        Material::Silicon => 1,
        Material::SiliconDioxide => 2,
        Material::SiliconCarbide => 3,
        Material::GalliumNitride => 4,
        Material::Photoresist => 5,
        Material::Metal => 6,
        Material::Dopant(b) => DOPANT_ID_BASE + (b.min(8) as i32),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_children_are_contiguous_for_shader_traversal() {
        let mut wafer = WaferState::new_substrate(300.0, 500.0, 10_000_000.0);
        wafer.set_region(
            Bounds { min: wafer.bounds.min, max: wafer.bounds.center() },
            Material::SiliconDioxide,
        );
        let flat = flatten_for_gpu(&wafer);
        assert_eq!(flat[0].material_id, -1);
        let first = flat[0].first_child as usize;
        assert!(first > 0);
        assert!(first + 8 <= flat.len());
        // Every direct child occupies the contiguous block the WGSL shader indexes.
        assert_eq!(flat[first + 0].min[0], wafer.bounds.min.0 as f32);
        assert_eq!(flat[first + 1].max[0], wafer.bounds.max.0 as f32);
    }

    #[test]
    fn gpu_record_matches_wgsl_storage_stride() {
        assert_eq!(std::mem::size_of::<GpuSvoNode>(), 48);
        assert_eq!(std::mem::align_of::<GpuSvoNode>(), 16);
    }

    #[test]
    fn single_leaf_wafer_flattens_to_one_node() {
        let wafer = WaferState::new_substrate(300.0, 500.0, 1.0);
        let flat = flatten_for_gpu(&wafer);
        assert_eq!(flat.len(), 1);
        assert_eq!(flat[0].material_id, material_to_id(Material::Silicon));
        assert_eq!(flat[0].first_child, -1);
    }
}
