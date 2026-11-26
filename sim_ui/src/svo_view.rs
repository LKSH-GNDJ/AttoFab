//! GPU presentation of the wafer: a fullscreen raymarch of the flattened SVO tile
//! plus an analytic bulk. All camera / palette decisions arrive pre-packed from
//! [`crate::viewport_math::pack_uniform`] so the CPU and the shader agree.

use bytemuck::{Pod, Zeroable};
use eframe::{egui, egui_wgpu::{self, wgpu}};
use eframe::egui_wgpu::wgpu::util::DeviceExt;
use sim_core::render_wgpu::GpuSvoNode;
use std::sync::Arc;

use crate::viewport_math::UNIFORM_VEC4S;

pub const SHADER_SOURCE: &str = include_str!("../../assets/shaders/wafer_render.wgsl");

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct CameraUniform { data: [[f32; 4]; UNIFORM_VEC4S] }

pub struct SvoRenderResources {
    pipeline: wgpu::RenderPipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    camera_buffer: wgpu::Buffer,
    node_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    uploaded_revision: u64,
}

impl SvoRenderResources {
    fn update(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, nodes: &[GpuSvoNode], revision: u64, camera: CameraUniform) {
        if revision != self.uploaded_revision {
            let bytes = bytemuck::cast_slice(nodes);
            if self.node_buffer.size() < bytes.len() as u64 {
                let size = (bytes.len() as u64).next_power_of_two().max(48);
                self.node_buffer = device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("attofab-svo-nodes"), size, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
                });
                self.bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("attofab-svo-bind-group"), layout: &self.bind_group_layout,
                    entries: &[
                        wgpu::BindGroupEntry { binding: 0, resource: self.camera_buffer.as_entire_binding() },
                        wgpu::BindGroupEntry { binding: 1, resource: self.node_buffer.as_entire_binding() },
                    ],
                });
            }
            if !bytes.is_empty() { queue.write_buffer(&self.node_buffer, 0, bytes); }
            self.uploaded_revision = revision;
        }
        queue.write_buffer(&self.camera_buffer, 0, bytemuck::bytes_of(&camera));
    }

    fn paint<'a>(&'a self, pass: &mut wgpu::RenderPass<'a>) {
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, &self.bind_group, &[]);
        pass.draw(0..3, 0..1);
    }
}

/// Build the pipeline and buffers (shared by the live viewport and the offscreen renderer).
pub fn create_resources(device: &wgpu::Device, format: wgpu::TextureFormat) -> SvoRenderResources {
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("attofab-svo-raycaster"), source: wgpu::ShaderSource::Wgsl(SHADER_SOURCE.into()) });
    let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("attofab-svo-bind-layout"),
        entries: &[
            wgpu::BindGroupLayoutEntry { binding: 0, visibility: wgpu::ShaderStages::FRAGMENT, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None }, count: None },
            wgpu::BindGroupLayoutEntry { binding: 1, visibility: wgpu::ShaderStages::FRAGMENT, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: true }, has_dynamic_offset: false, min_binding_size: None }, count: None },
        ],
    });
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("attofab-svo-pipeline-layout"), bind_group_layouts: &[&layout], push_constant_ranges: &[] });
    let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("attofab-svo-pipeline"), layout: Some(&pipeline_layout),
        vertex: wgpu::VertexState { module: &shader, entry_point: "vs_main", buffers: &[] },
        fragment: Some(wgpu::FragmentState { module: &shader, entry_point: "fs_main", targets: &[Some(wgpu::ColorTargetState { format, blend: Some(wgpu::BlendState::ALPHA_BLENDING), write_mask: wgpu::ColorWrites::ALL })] }),
        primitive: wgpu::PrimitiveState::default(), depth_stencil: None, multisample: wgpu::MultisampleState::default(), multiview: None,
    });
    let camera = CameraUniform { data: [[0.0; 4]; UNIFORM_VEC4S] };
    let camera_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("attofab-camera"), contents: bytemuck::bytes_of(&camera), usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST });
    let node_buffer = device.create_buffer(&wgpu::BufferDescriptor { label: Some("attofab-svo-nodes"), size: 48, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("attofab-svo-bind-group"), layout: &layout,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: camera_buffer.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: node_buffer.as_entire_binding() },
        ],
    });
    SvoRenderResources { pipeline, bind_group_layout: layout, camera_buffer, node_buffer, bind_group, uploaded_revision: u64::MAX }
}

/// Render one frame with the real shader into an RGBA8 image (headless / regression use).
pub fn render_offscreen(device: &wgpu::Device, queue: &wgpu::Queue, nodes: &[GpuSvoNode], uniform: [[f32; 4]; UNIFORM_VEC4S], w: u32, h: u32, background: [f64; 4]) -> Vec<u8> {
    render_offscreen_in(device, queue, nodes, uniform, (w, h), (0, 0, w, h), background)
}

/// Like [`render_offscreen`] but draws into a larger target through an *inset* viewport
/// `(x, y, w, h)`, exactly like the live window does. Returns the whole target.
pub fn render_offscreen_in(device: &wgpu::Device, queue: &wgpu::Queue, nodes: &[GpuSvoNode], mut uniform: [[f32; 4]; UNIFORM_VEC4S], target: (u32, u32), vp: (u32, u32, u32, u32), background: [f64; 4]) -> Vec<u8> {
    let (tw, th) = target;
    uniform[4][0] = vp.2 as f32; uniform[4][1] = vp.3 as f32; uniform[4][3] = vp.0 as f32; uniform[5][3] = vp.1 as f32;
    let format = wgpu::TextureFormat::Rgba8Unorm;
    let mut res = create_resources(device, format);
    res.update(device, queue, nodes, 1, CameraUniform { data: uniform });
    let tex = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("offscreen"), size: wgpu::Extent3d { width: tw, height: th, depth_or_array_layers: 1 }, mip_level_count: 1, sample_count: 1,
        dimension: wgpu::TextureDimension::D2, format, usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC, view_formats: &[],
    });
    let view = tex.create_view(&wgpu::TextureViewDescriptor::default());
    let row = (tw * 4 + 255) / 256 * 256;
    let out = device.create_buffer(&wgpu::BufferDescriptor { label: Some("readback"), size: (row * th) as u64, usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ, mapped_at_creation: false });
    let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("offscreen") });
    {
        let mut pass = enc.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("offscreen-pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment { view: &view, resolve_target: None, ops: wgpu::Operations { load: wgpu::LoadOp::Clear(wgpu::Color { r: background[0], g: background[1], b: background[2], a: background[3] }), store: wgpu::StoreOp::Store } })],
            depth_stencil_attachment: None, timestamp_writes: None, occlusion_query_set: None,
        });
        pass.set_viewport(vp.0 as f32, vp.1 as f32, vp.2 as f32, vp.3 as f32, 0.0, 1.0);
        pass.set_scissor_rect(vp.0, vp.1, vp.2, vp.3);
        res.paint(&mut pass);
    }
    enc.copy_texture_to_buffer(wgpu::ImageCopyTexture { texture: &tex, mip_level: 0, origin: wgpu::Origin3d::ZERO, aspect: wgpu::TextureAspect::All },
        wgpu::ImageCopyBuffer { buffer: &out, layout: wgpu::ImageDataLayout { offset: 0, bytes_per_row: Some(row), rows_per_image: Some(th) } }, wgpu::Extent3d { width: tw, height: th, depth_or_array_layers: 1 });
    queue.submit(Some(enc.finish()));
    let slice = out.slice(..);
    slice.map_async(wgpu::MapMode::Read, |_| {});
    device.poll(wgpu::Maintain::Wait);
    let data = slice.get_mapped_range();
    let mut img = Vec::with_capacity((tw * th * 4) as usize);
    for y in 0..th as usize { img.extend_from_slice(&data[y * row as usize..y * row as usize + (tw * 4) as usize]); }
    img
}

#[derive(Clone)]
pub struct SvoViewportCallback {
    /// Tile-relative flattened SVO. Re-uploaded only when `revision` changes.
    pub nodes: Arc<Vec<GpuSvoNode>>,
    pub revision: u64,
    pub uniform: [[f32; 4]; UNIFORM_VEC4S],
    /// Viewport rectangle in egui points; converted to framebuffer pixels in `prepare`.
    pub rect: egui::Rect,
}

impl SvoViewportCallback {
    pub fn install(cc: &eframe::CreationContext<'_>) -> bool {
        let Some(state) = cc.wgpu_render_state.as_ref() else { return false; };
        let res = create_resources(&state.device, state.target_format);
        state.renderer.write().callback_resources.insert(res);
        true
    }

    pub fn add_to(&self, ui: &mut egui::Ui, rect: egui::Rect) {
        ui.painter().add(egui_wgpu::Callback::new_paint_callback(rect, self.clone()));
    }
}

impl egui_wgpu::CallbackTrait for SvoViewportCallback {
    fn prepare(&self, device: &wgpu::Device, queue: &wgpu::Queue, screen: &egui_wgpu::ScreenDescriptor, _encoder: &mut wgpu::CommandEncoder, resources: &mut egui_wgpu::CallbackResources) -> Vec<wgpu::CommandBuffer> {
        // The shader works in framebuffer pixels: give it the viewport's real origin and size.
        let (origin, size) = crate::viewport_math::viewport_pixels([self.rect.min.x, self.rect.min.y], [self.rect.width(), self.rect.height()], screen.pixels_per_point);
        let mut u = self.uniform;
        u[4][0] = size[0]; u[4][1] = size[1]; u[4][3] = origin[0]; u[5][3] = origin[1];
        if let Some(r) = resources.get_mut::<SvoRenderResources>() { r.update(device, queue, &self.nodes, self.revision, CameraUniform { data: u }); }
        Vec::new()
    }

    fn paint<'a>(&'a self, _info: egui::PaintCallbackInfo, pass: &mut wgpu::RenderPass<'a>, resources: &'a egui_wgpu::CallbackResources) {
        if let Some(r) = resources.get::<SvoRenderResources>() { r.paint(pass); }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The shader must parse and validate, and its uniform block must be exactly the
    /// 48 × vec4 the CPU packs (768 bytes).
    #[test]
    fn wgsl_parses_validates_and_matches_the_uniform_layout() {
        let module = naga::front::wgsl::parse_str(SHADER_SOURCE).expect("WGSL parse");
        let info = naga::valid::Validator::new(naga::valid::ValidationFlags::all(), naga::valid::Capabilities::all()).validate(&module).expect("WGSL validation");
        let _ = info;
        let (_, ty) = module.types.iter().find(|(_, t)| t.name.as_deref() == Some("CameraUniform")).expect("CameraUniform struct");
        match &ty.inner {
            naga::TypeInner::Struct { span, .. } => assert_eq!(*span as usize, UNIFORM_VEC4S * 16),
            other => panic!("unexpected {other:?}"),
        }
        assert_eq!(std::mem::size_of::<CameraUniform>(), UNIFORM_VEC4S * 16);
        assert!(module.entry_points.iter().any(|e| e.name == "vs_main") && module.entry_points.iter().any(|e| e.name == "fs_main"));
    }

    #[test]
    fn gpu_node_layout_matches_the_shader_struct() {
        // SvoNode { bmin: vec3, bmax: vec3, material_id: i32, first_child: i32 } is 48 bytes in a storage array.
        assert_eq!(std::mem::size_of::<GpuSvoNode>(), 48);
    }
}
