//! `gds_info file.gds` — what is in this layout? Cells, hierarchy, units, layers, feature sizes.
use sim_core::io::gds::{parse, summarize};

fn main() {
    let Some(path) = std::env::args().nth(1) else { eprintln!("usage: gds_info file.gds"); std::process::exit(2) };
    let bytes = std::fs::read(&path).unwrap_or_else(|e| { eprintln!("{path}: {e}"); std::process::exit(2) });
    let lay = match parse(&bytes) { Ok(l) => l, Err(e) => { eprintln!("not loadable: {e}"); std::process::exit(1) } };
    println!("file: {path}\nbytes: {}\ndatabase unit: {} nm\ncells: {}", bytes.len(), lay.db_unit_m * 1e9, lay.structures.len());
    println!("top cells: {:?} (default: {:?})", lay.top_cells(), lay.default_top());
    let flat = lay.flatten(None).unwrap_or_else(|e| { eprintln!("flatten: {e}"); std::process::exit(1) });
    println!("flattened polygons: {}", flat.len());
    println!("{:>6} {:>7} {:>10} {:>26} {:>16} {:>16}", "layer", "dtype", "polygons", "bbox size (µm)", "p10 feature nm", "median feature nm");
    for l in summarize(&flat) {
        let (w, h) = l.size_nm();
        println!("{:>6} {:>7} {:>10} {:>12.1} x {:<11.1} {:>16.0} {:>16.0}", l.layer, l.datatype, l.polygons, w / 1e3, h / 1e3, l.feature_p10_nm, l.feature_med_nm);
    }
}
