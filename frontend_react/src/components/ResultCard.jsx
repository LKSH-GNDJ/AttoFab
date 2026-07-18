export default function ResultCard({ result }) {
  if (!result) return null;
  const { run_id, summary } = result;

  return (
    <div className="result-card">
      <h3>Run #{run_id}</h3>

      <div className="result-section">
        <h4>Material composition</h4>
        <table>
          <tbody>
            {Object.entries(summary.material_pct).map(([mat, pct]) => (
              <tr key={mat}>
                <td>{mat}</td>
                <td>{pct}%</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <div className="result-section">
        <h4>Oxide thickness</h4>
        <table>
          <tbody>
            <tr><td>Average</td><td>{summary.oxide_nm.avg_nm} nm</td></tr>
            <tr><td>Max</td><td>{summary.oxide_nm.max_nm} nm</td></tr>
            <tr><td>Min</td><td>{summary.oxide_nm.min_nm} nm</td></tr>
          </tbody>
        </table>
      </div>

      {Object.keys(summary.species_peak_cm3 || {}).length > 0 && (
        <div className="result-section">
          <h4>Peak dopant concentration</h4>
          <table>
            <tbody>
              {Object.entries(summary.species_peak_cm3).map(([dopant, peak]) => (
                <tr key={dopant}>
                  <td>{dopant}</td>
                  <td>{peak} /cm&sup3;</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      <div className="result-section">
        <h4>Process log</h4>
        <ul className="process-log">
          {summary.process_steps.map((step, i) => (
            <li key={i}>{step}</li>
          ))}
        </ul>
      </div>
    </div>
  );
}
