//! AttoFab native Android platform shell.
//!
//! Built as a `cdylib` (`libapp_android.so`) and loaded by Android's
//! `NativeActivity` via the `android_activity` crate's glue, which calls
//! `android_main` below with an `AndroidApp` handle wrapping the
//! `ANativeActivity`/`ANativeWindow` the app has been given.
//!
//! As with `app_web`, all session/state logic lives in the
//! platform-agnostic `bootstrap` module so it stays testable on a plain
//! host target; this file is just the NDK-facing glue.

pub mod bootstrap;

#[cfg(target_os = "android")]
mod android_entry {
    use super::bootstrap::Session;
    use android_activity::{AndroidApp, MainEvent, PollEvent};

    #[no_mangle]
    fn android_main(app: AndroidApp) {
        android_logger::init_once(
            android_logger::Config::default().with_max_level(log::LevelFilter::Info),
        );
        log::info!("AttoFab (android) starting — GNU GPL-3.0-or-later");

        let mut session = Session::new_default();
        log::info!("{}", session.status_line());

        // A full graphical build creates a `wgpu::Surface` from the
        // `ANativeWindow` handed over in `PollEvent::Main(MainEvent::InitWindow
        // { .. })` below, then drives the same `sim_ui` panel logic used by
        // the desktop/web shells against it (touch input replacing
        // mouse/keyboard for pan/zoom/parameter edits). That wgpu wiring is
        // left as documentation here, matching the same scope boundary as
        // `app_web`'s entry point.
        let mut running = true;
        let mut window_ready = false;

        while running {
            app.poll_events(Some(std::time::Duration::from_millis(16)), |event| {
                match event {
                    PollEvent::Main(MainEvent::InitWindow { .. }) => {
                        window_ready = true;
                        log::info!("Native window ready; ripe for a wgpu surface here.");
                    }
                    PollEvent::Main(MainEvent::TerminateWindow { .. }) => {
                        window_ready = false;
                    }
                    PollEvent::Main(MainEvent::Destroy) => {
                        running = false;
                    }
                    _ => {}
                }
            });

            if window_ready {
                // Placeholder per-frame tick: in a full build this steps
                // `session.engine.poll_results()` and redraws the wgpu
                // surface. Marking a stage dirty here is a stand-in for a
                // real touch-driven parameter-table edit event.
                let _ = &mut session;
            }
        }
    }
}
