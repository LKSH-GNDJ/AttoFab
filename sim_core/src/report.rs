//! Process-traceability report (PRD "complete process traceability"): the recipe,
//! every stage's parameters, input/output wafer fingerprints, fidelity,
//! sub-steps, metrics, warnings and the metrology summary, as Markdown.

use crate::measure::Status;
use crate::params::Severity;
use crate::session::{Session, StageStatus};
use std::fmt::Write;

pub fn markdown_report(s: &Session) -> String {
    let mut o = String::new();
    let st = &s.state;
    let w = &st.recipe.window;
    let _ = writeln!(o, "# AttoFab process report\n");
    let _ = writeln!(o, "- Wafer: Ø{} mm × {} µm, octree cell {} nm, high-resolution tile {} nm", st.spec.diameter_mm, st.spec.thickness_um, st.spec.leaf_nm, st.spec.tile_edge_nm());
    let _ = writeln!(o, "- Process site (wafer-local): ({:.0}, {:.0}) nm, window {:.0} × {:.0} nm", w.center[0], w.center[1], w.size[0], w.size[1]);
    let _ = writeln!(o, "- Seed: {} · committed revision: {} · initial wafer fingerprint `{:016x}`\n", st.recipe.seed, st.revision, st.initial.fingerprint());
    let statuses = s.statuses();
    for (i, rs) in st.recipe.stages.iter().enumerate() {
        let status = &statuses[i].1;
        let _ = writeln!(o, "## {}. {} — {}\n", i + 1, rs.params.stage().label(), match status { StageStatus::Committed => "committed".to_string(), StageStatus::Stale => "STALE (an upstream change invalidated this result)".into(), StageStatus::NotRun => "not run".into(), StageStatus::Disabled => "disabled".into(), StageStatus::Failed(m) => format!("FAILED: {m}") });
        let _ = writeln!(o, "Parameters: `{:?}`\n", rs.params);
        if let (StageStatus::Committed, Some(r)) = (status, st.records.get(&rs.id)) {
            let _ = writeln!(o, "- rev {} · fidelity **{}** · input `{:016x}` → output `{:016x}` · nominal duration {:.0} s", r.revision, r.outcome.fidelity.label(), r.input_fp, r.output_fp, r.outcome.duration_s);
            let _ = writeln!(o, "- geometry changed: {}", r.outcome.geometry_changed);
            let _ = writeln!(o, "- sub-steps: {}", r.outcome.substeps.iter().map(|x| format!("{} ({})", x.name, x.kind.label())).collect::<Vec<_>>().join(" → "));
            for m in &r.outcome.metrics { let _ = writeln!(o, "- {}: {:.4} {}", m.name, m.value, m.unit); }
            for i in &r.outcome.issues { let _ = writeln!(o, "- {}: {}", if i.severity == Severity::Error { "ERROR" } else { "warning" }, i.message); }
            if let Some(rep) = &r.outcome.report {
                let _ = writeln!(o, "\n| Measurement | Value | Nominal ± tol | Status |\n|---|---|---|---|");
                for m in &rep.measurements {
                    let _ = writeln!(o, "| {} | {} | {:.1} ± {:.1} {} | {} |", m.name, m.value.map(|v| format!("{v:.2}")).unwrap_or("—".into()), m.nominal, m.tolerance, m.unit, match m.status { Status::Pass => "pass", Status::Fail => "FAIL", Status::NotMeasured => "not measured" });
                }
                let _ = writeln!(o, "\nMeasured wafer fingerprint `{:016x}` · recipe: **{}**", rep.wafer_fingerprint, match rep.pass { Some(true) => "PASS", Some(false) => "FAIL", None => "no measurable targets" });
            }
        }
        let _ = writeln!(o);
    }
    o
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn report_names_every_stage_fidelity_and_the_measured_fingerprint() {
        let mut s = Session::with_default_recipe();
        s.run_recipe();
        let md = markdown_report(&s);
        for l in ["Substrate Prep", "Deposition", "Lithography", "Etching", "Doping & Anneal", "Final Clean", "Metrology"] { assert!(md.contains(l), "{l}"); }
        assert!(md.contains("Approximated") && md.contains("recipe: **PASS**"));
        let fp = s.committed_wafer().fingerprint();
        assert!(md.contains(&format!("{fp:016x}")), "report must cite the final wafer fingerprint");
        let mut s2 = Session::with_default_recipe();
        assert!(markdown_report(&s2).contains("not run"));
        s2.run_recipe();
    }
}
