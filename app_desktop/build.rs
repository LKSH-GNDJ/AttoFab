// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Embeds the AttoFab icon (and version info) into the Windows executable so Explorer, the taskbar and
//! Alt-Tab show it. Never fails the build: without the Windows resource compiler it only prints a warning.
fn main() {
    println!("cargo:rerun-if-changed=assets/attofab.ico");
    #[cfg(windows)]
    {
        let mut res = winresource::WindowsResource::new();
        res.set_icon("assets/attofab.ico");
        res.set("ProductName", "AttoFab");
        res.set("FileDescription", "AttoFab — Semiconductor Fabrication Simulator");
        if let Err(e) = res.compile() { println!("cargo:warning=could not embed the exe icon: {e}"); }
    }
}
