// AttoFab MK1 SVO raycaster (tile + implicit bulk).
//
// Coordinates are RELATIVE to the tile origin (the centre of the high-
// resolution tile at the substrate surface), computed in f64 on the CPU before
// the f32 cast. That keeps nanometre features exact even when the tile sits
// 100 mm from the wafer centre.
//
// Scene = flattened SVO tile + an analytic bulk (a silicon cylinder that is
// subtracted inside the tile box). The camera basis, field of view and the
// colour palette all come from the CPU so picking, overlays and rendering share
// one definition.

struct SvoNode {
    bmin: vec3<f32>,
    bmax: vec3<f32>,
    material_id: i32,
    first_child: i32,
};

// 48 x vec4 = 768 bytes (34 used + padding so the Rust array has a bytemuck-supported length);
// mirrors Rust `CameraUniform { data: [[f32; 4]; 48] }`.
struct CameraUniform {
    camera_pos: vec4<f32>,       //  0
    camera_forward: vec4<f32>,   //  1
    camera_right: vec4<f32>,     //  2
    camera_up: vec4<f32>,        //  3
    viewport: vec4<f32>,         //  4  (width px, height px, tan(half fov), origin x px)
    view_slice: vec4<f32>,       //  5  (view mode, slice z, clip enabled, origin y px)
    clip_plane: vec4<f32>,       //  6  keep dot(n, p) + d >= 0
    domain: vec4<f32>,           //  7  (radius, centre x, centre y, bulk enabled)
    bulk_z: vec4<f32>,           //  8  (bottom, top, 0, 0)
    tile_min: vec4<f32>,         //  9
    tile_max: vec4<f32>,         // 10
    selected_min: vec4<f32>,     // 11
    selected_max: vec4<f32>,     // 12
    selection_enabled: vec4<f32>,// 13
    visible: array<vec4<f32>, 4>,    // 14..17
    palette: array<vec4<f32>, 16>,   // 18..33
    reserved: array<vec4<f32>, 14>,  // 34..47
};

@group(0) @binding(0) var<uniform> camera: CameraUniform;
@group(0) @binding(1) var<storage, read> svo_nodes: array<SvoNode>;

const NONE: vec2<f32> = vec2<f32>(1e30, -1e30);

fn slot_visible(id: i32) -> bool {
    let v = camera.visible[id / 4];
    let k = id % 4;
    if (k == 0) { return v.x > 0.5; }
    if (k == 1) { return v.y > 0.5; }
    if (k == 2) { return v.z > 0.5; }
    return v.w > 0.5;
}

fn ray_box(ro: vec3<f32>, rd_inv: vec3<f32>, bmin: vec3<f32>, bmax: vec3<f32>) -> vec2<f32> {
    let t0 = (bmin - ro) * rd_inv;
    let t1 = (bmax - ro) * rd_inv;
    let lo = min(t0, t1);
    let hi = max(t0, t1);
    return vec2<f32>(max(max(lo.x, lo.y), lo.z), min(min(hi.x, hi.y), hi.z));
}

// Interval inside the wafer's circular XY domain.
fn ray_disk(ro: vec3<f32>, rd: vec3<f32>) -> vec2<f32> {
    let radius = max(camera.domain.x, 1.0);
    let p = ro.xy - camera.domain.yz;
    let a = dot(rd.xy, rd.xy);
    let c = dot(p, p) - radius * radius;
    if (a < 1e-10) {
        if (c <= 0.0) { return vec2<f32>(-1e30, 1e30); }
        return NONE;
    }
    let b = 2.0 * dot(p, rd.xy);
    let disc = b * b - 4.0 * a * c;
    if (disc < 0.0) { return NONE; }
    let r = sqrt(disc);
    return vec2<f32>((-b - r) / (2.0 * a), (-b + r) / (2.0 * a));
}

// Interval inside the bulk's z slab.
fn ray_slab(ro: vec3<f32>, rd: vec3<f32>) -> vec2<f32> {
    let lo = camera.bulk_z.x;
    let hi = camera.bulk_z.y;
    if (abs(rd.z) < 1e-9) {
        if (ro.z >= lo && ro.z <= hi) { return vec2<f32>(-1e30, 1e30); }
        return NONE;
    }
    let a = (lo - ro.z) / rd.z;
    let b = (hi - ro.z) / rd.z;
    return vec2<f32>(min(a, b), max(a, b));
}

// Section plane (cross-section views): keep dot(n, p) + d >= 0.
fn clip_interval(ro: vec3<f32>, rd: vec3<f32>, t: vec2<f32>) -> vec2<f32> {
    if (camera.view_slice.z < 0.5) { return t; }
    let n = camera.clip_plane.xyz;
    let num = dot(n, ro) + camera.clip_plane.w;
    let den = dot(n, rd);
    if (abs(den) < 1e-9) {
        if (num < 0.0) { return NONE; }
        return t;
    }
    let tc = -num / den;
    if (den > 0.0) { return vec2<f32>(max(t.x, tc), t.y); }
    return vec2<f32>(t.x, min(t.y, tc));
}

fn blend(accum: vec4<f32>, c: vec4<f32>) -> vec4<f32> {
    return vec4<f32>(accum.rgb + (1.0 - accum.a) * c.a * c.rgb, accum.a + (1.0 - accum.a) * c.a);
}

const MAX_STEPS: i32 = 6000;
const STACK_DEPTH: i32 = 256;

fn raymarch_tile(ro: vec3<f32>, rd: vec3<f32>, window: vec2<f32>) -> vec4<f32> {
    let rd_inv = 1.0 / rd;
    var mask: u32 = 0u;
    if (rd.x < 0.0) { mask = mask | 1u; }
    if (rd.y < 0.0) { mask = mask | 2u; }
    if (rd.z < 0.0) { mask = mask | 4u; }
    var accum = vec4<f32>(0.0);
    var stack: array<i32, STACK_DEPTH>;
    var sp: i32 = 1;
    stack[0] = 0;
    var steps: i32 = 0;
    loop {
        if (sp <= 0 || steps >= MAX_STEPS || accum.a >= 0.98) { break; }
        steps = steps + 1;
        sp = sp - 1;
        let node = svo_nodes[stack[sp]];
        if (camera.view_slice.x > 3.5 && (node.bmax.z < camera.view_slice.y || node.bmin.z > camera.view_slice.y)) { continue; }
        let box_hit = ray_box(ro, rd_inv, node.bmin, node.bmax);
        let hit = vec2<f32>(max(box_hit.x, window.x), min(box_hit.y, window.y));
        if (hit.x > hit.y || hit.y < 0.0) { continue; }
        if (node.material_id >= 0) {
            if (node.material_id == 0 || !slot_visible(node.material_id)) { continue; }
            var c = camera.palette[node.material_id];
            if (camera.selection_enabled.x > 0.5) {
                let same_min = all(abs(node.bmin - camera.selected_min.xyz) < vec3<f32>(0.01));
                let same_max = all(abs(node.bmax - camera.selected_max.xyz) < vec3<f32>(0.01));
                if (same_min && same_max) { c = vec4<f32>(1.0, 0.82, 0.18, 0.9); }
            }
            accum = blend(accum, c);
        } else if (node.first_child >= 0 && sp + 8 <= STACK_DEPTH) {
            // Push far-to-near so the nearest octant is popped first.
            for (var k: i32 = 7; k >= 0; k = k - 1) {
                stack[sp] = node.first_child + i32(u32(k) ^ mask);
                sp = sp + 1;
            }
        }
    }
    return accum;
}

fn trace(ro: vec3<f32>, rd: vec3<f32>) -> vec4<f32> {
    let rd_inv = 1.0 / rd;
    var window = ray_disk(ro, rd);
    window = clip_interval(ro, rd, window);

    // Tile box interval (also limited by the wafer disk and the section plane).
    let tb = ray_box(ro, rd_inv, camera.tile_min.xyz, camera.tile_max.xyz);
    let tile = vec2<f32>(max(tb.x, window.x), min(tb.y, window.y));
    let tile_hit = tile.x <= tile.y && tile.y >= 0.0;

    // Bulk interval: disk ∩ slab ∩ section. In slice mode it is cut at the slice height.
    var bulk = ray_slab(ro, rd);
    bulk = vec2<f32>(max(bulk.x, window.x), min(bulk.y, window.y));
    var bulk_on = camera.domain.w > 0.5 && bulk.x <= bulk.y && bulk.y >= 0.0 && slot_visible(1);
    if (camera.view_slice.x > 3.5 && (camera.view_slice.y < camera.bulk_z.x || camera.view_slice.y > camera.bulk_z.y)) { bulk_on = false; }
    let bulk_col = camera.palette[1];

    var pre_end = bulk.y;
    var post_start = bulk.y + 1.0;
    if (tile_hit) {
        pre_end = min(bulk.y, tb.x);
        post_start = max(bulk.x, tb.y);
    }
    // Bulk in front of the tile is opaque silicon.
    if (bulk_on && bulk.x < pre_end && pre_end >= 0.0) { return vec4<f32>(bulk_col.rgb, 1.0); }

    var accum = vec4<f32>(0.0);
    if (tile_hit) { accum = raymarch_tile(ro, rd, window); }
    if (bulk_on && tile_hit && post_start < bulk.y && accum.a < 0.98) {
        accum = blend(accum, vec4<f32>(bulk_col.rgb, 1.0));
    }
    return accum;
}

@vertex
fn vs_main(@builtin(vertex_index) i: u32) -> @builtin(position) vec4<f32> {
    var p = array<vec2<f32>, 3>(vec2<f32>(-1.0, -1.0), vec2<f32>(3.0, -1.0), vec2<f32>(-1.0, 3.0));
    return vec4<f32>(p[i], 0.0, 1.0);
}

@fragment
fn fs_main(@builtin(position) frag: vec4<f32>) -> @location(0) vec4<f32> {
    // @builtin(position) is in FRAMEBUFFER pixels; the viewport sits inside the window, so subtract its origin.
    let uv = (frag.xy - vec2<f32>(camera.viewport.w, camera.view_slice.w)) / camera.viewport.xy;
    let ndc = vec2<f32>(uv.x * 2.0 - 1.0, 1.0 - uv.y * 2.0);
    let aspect = camera.viewport.x / max(camera.viewport.y, 1.0);
    let t = camera.viewport.z;
    let rd = normalize(camera.camera_forward.xyz + camera.camera_right.xyz * ndc.x * aspect * t + camera.camera_up.xyz * ndc.y * t);
    let c = trace(camera.camera_pos.xyz, rd);
    if (c.a < 0.01) { discard; }
    return c;
}
