// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Solve the 3-D particle-flux field of a trench or a hole and compare the bottom coverage
//! with the two-parameter step-coverage model (PRD roadmap 6).
//!
//! ```text
//! cargo run --release -p atto_cli --bin field3d -- --shape hole --width 16 --depth 32 --sticking 0.3
//! ```
//!
//! Options (all optional): `--shape trench|hole`, `--width N` and `--depth N` in cells,
//! `--cell-nm X`, `--sticking S`, `--rays N`, `--threads N`, `--seed N`, and for a tilted beam
//! instead of a cosine source `--tilt-deg X` with `--azimuth-deg X` and `--spread-deg X`.

use atto_process::deposition::Precursor;
use atto_process::field3d::{Angular, Config, Partial, Plan, Source3d, VoxelGrid};
use atto_units::Nm;

struct Args {
    hole: bool,
    width: usize,
    depth: usize,
    cell_nm: f64,
    sticking: f64,
    rays: u64,
    threads: usize,
    seed: u64,
    tilt_deg: Option<f64>,
    azimuth_deg: f64,
    spread_deg: f64,
}

fn number<T: std::str::FromStr>(flag: &str, value: Option<&String>) -> Result<T, String> {
    let v = value.ok_or_else(|| format!("{flag} needs a value"))?;
    v.parse::<T>().map_err(|_| format!("{flag}: '{v}' is not a valid number"))
}

fn parse(argv: &[String]) -> Result<Args, String> {
    let mut a = Args {
        hole: false,
        width: 16,
        depth: 32,
        cell_nm: 10.0,
        sticking: 1.0,
        rays: 400_000,
        threads: 1,
        seed: 1,
        tilt_deg: None,
        azimuth_deg: 0.0,
        spread_deg: 0.0,
    };
    let mut i = 0;
    while i < argv.len() {
        let flag = argv[i].as_str();
        let value = argv.get(i + 1);
        match flag {
            "--shape" => {
                a.hole = match value.map(|s| s.as_str()) {
                    Some("hole") => true,
                    Some("trench") => false,
                    _ => return Err("--shape must be 'trench' or 'hole'".into()),
                }
            }
            "--width" => a.width = number(flag, value)?,
            "--depth" => a.depth = number(flag, value)?,
            "--cell-nm" => a.cell_nm = number(flag, value)?,
            "--sticking" => a.sticking = number(flag, value)?,
            "--rays" => a.rays = number(flag, value)?,
            "--threads" => a.threads = number(flag, value)?,
            "--seed" => a.seed = number(flag, value)?,
            "--tilt-deg" => a.tilt_deg = Some(number(flag, value)?),
            "--azimuth-deg" => a.azimuth_deg = number(flag, value)?,
            "--spread-deg" => a.spread_deg = number(flag, value)?,
            other => return Err(format!("unknown option '{other}'")),
        }
        i += 2;
    }
    Ok(a)
}

fn run(argv: &[String]) -> Result<(), String> {
    let a = parse(argv)?;
    if a.width == 0 || a.depth == 0 {
        return Err("--width and --depth must be at least 1 cell".into());
    }
    let floor = 4usize;
    let nx = (a.width * 3).max(24);
    let ny = if a.hole { nx } else { 4 };
    let nz = floor + a.depth + 8;
    let cell = Nm::new(a.cell_nm);
    let built = if a.hole {
        VoxelGrid::hole(nx, ny, nz, cell, a.width, a.depth, floor)
    } else {
        VoxelGrid::trench(nx, ny, nz, cell, a.width, a.depth, floor)
    };
    let grid = built.map_err(|e| e.to_string())?;

    let source = match a.tilt_deg {
        Some(t) => Source3d { angular: Angular::Beam, tilt_deg: t, azimuth_deg: a.azimuth_deg, spread_deg: a.spread_deg },
        None => Source3d::default(),
    };
    let cfg = Config { rays: a.rays, seed: a.seed, sticking: a.sticking, source };
    if a.threads == 0 || a.threads > 64 {
        return Err("--threads must be between 1 and 64".into());
    }
    let plan = Plan::new(&grid, &cfg).map_err(|e| e.to_string())?;
    let chunks = plan.chunk_count();
    let workers = a.threads;
    let mut parts: Vec<Partial> = Vec::new();
    if workers == 1 {
        for c in 0..chunks {
            parts.push(plan.run_chunk(c));
        }
    } else {
        // the library starts no threads; this binary does. Chunks are independent and
        // `finish` orders them, so the result is the same as the single-threaded run.
        std::thread::scope(|scope| {
            let plan = &plan;
            let handles: Vec<_> = (0..workers)
                .map(|w| {
                    scope.spawn(move || {
                        let mut mine = Vec::new();
                        let mut c = w;
                        while c < chunks {
                            mine.push(plan.run_chunk(c));
                            c += workers;
                        }
                        mine
                    })
                })
                .collect();
            for h in handles {
                parts.extend(h.join().expect("a flux worker panicked"));
            }
        });
    }
    let sol = plan.finish(parts).map_err(|e| e.to_string())?;

    let xs = (nx / 2 - 1)..(nx / 2 + 1);
    let ys = if a.hole { (ny / 2 - 1)..(ny / 2 + 1) } else { 0..ny };
    let field = sol.mean_top_flux(0..4, 0..4).ok_or("no field columns to average")?;
    let bottom = sol.mean_top_flux(xs, ys).ok_or("no bottom columns to average")?;
    let b = sol.budget();
    let aspect = a.depth as f64 / a.width as f64;

    println!("AttoFab 3-D flux solve");
    println!("  geometry      {} {} x {} cells ({} nm cells), aspect ratio {:.2}", if a.hole { "square hole" } else { "trench" }, a.width, a.depth, a.cell_nm, aspect);
    println!("  source        {}", match a.tilt_deg {
        Some(t) => format!("beam, tilt {t} deg toward {} deg, spread {} deg", a.azimuth_deg, a.spread_deg),
        None => "cosine-law (Lambertian)".to_string(),
    });
    println!("  rays          {} (seed {}, sticking {})", a.rays, a.seed, a.sticking);
    println!("  field flux    {field:.4}  (open surface, reference 1)");
    println!("  bottom flux   {bottom:.4}  (centre of the feature floor)");
    println!("  bottom/field  {:.4}  <- 3-D step coverage", bottom / field);
    match (Precursor { sticking: a.sticking, ar_half: 4.0 }).bottom_coverage(aspect) {
        Ok(c) => println!("  two-parameter model (ar_half 4.0) {c:.4}  <- the number the deposition stage uses today"),
        Err(e) => println!("  two-parameter model refused: {e}"),
    }
    println!("  weight budget launched {:.0}, absorbed {:.1}, escaped {:.1}, truncated {:.1} (imbalance {:.2e})", b.launched, b.absorbed, b.escaped, b.truncated, b.imbalance());
    if b.imbalance() > 1e-6 || b.truncated_fraction() > 0.01 {
        return Err("the weight budget does not balance or too much was truncated: do not trust this result".into());
    }
    Ok(())
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if let Err(e) = run(&argv) {
        eprintln!("field3d: {e}");
        std::process::exit(1);
    }
}
