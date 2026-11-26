//! Browser entry point. The same egui application is used by desktop and web.

#[cfg(target_arch = "wasm32")]
mod web_entry {
    use wasm_bindgen::prelude::*;

    #[wasm_bindgen(start)]
    pub fn start() -> Result<(), JsValue> {
        console_error_panic_hook::set_once();
        let _ = console_log::init_with_level(log::Level::Info);
        wasm_bindgen_futures::spawn_local(async move {
            let runner = eframe::WebRunner::new();
            let result = runner
                .start(
                    "attofab-canvas",
                    eframe::WebOptions::default(),
                    Box::new(|cc| Box::new(sim_ui::app::AttoFabApp::new(cc))),
                )
                .await;
            if let Err(err) = result {
                log::error!("Unable to start AttoFab web UI: {err}");
            }
        });
        Ok(())
    }
}

