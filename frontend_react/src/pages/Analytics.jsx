import { useEffect, useState } from 'react';
import { listRuns } from '../api/client.js';

export default function Analytics() {
  const [runs, setRuns] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    listRuns(200)
      .then(setRuns)
      .catch((e) => setError(e.message));
  }, []);

  if (error) return <div className="page"><p className="error-text">{error}</p></div>;
  if (!runs) return <div className="page"><p className="muted">Loading\u2026</p></div>;
  if (runs.length === 0) return <div className="page"><p className="muted">No runs yet.</p></div>;

  const avgOxide = (
    runs.reduce((sum, r) => sum + (r.summary.oxide_nm?.avg_nm || 0), 0) / runs.length
  ).toFixed(1);

  const materialTotals = {};
  for (const r of runs) {
    for (const [mat, pct] of Object.entries(r.summary.material_pct || {})) {
      materialTotals[mat] = (materialTotals[mat] || 0) + pct;
    }
  }
  const materialAvg = Object.fromEntries(
    Object.entries(materialTotals).map(([mat, total]) => [mat, (total / runs.length).toFixed(1)])
  );

  return (
    <div className="page">
      <h1>Analytics</h1>
      <div className="panel">
        <h3>Overview</h3>
        <table>
          <tbody>
            <tr><td>Total runs</td><td>{runs.length}</td></tr>
            <tr><td>Average oxide thickness</td><td>{avgOxide} nm</td></tr>
          </tbody>
        </table>
      </div>
      <div className="panel">
        <h3>Average material composition across all runs</h3>
        <table>
          <tbody>
            {Object.entries(materialAvg).map(([mat, pct]) => (
              <tr key={mat}><td>{mat}</td><td>{pct}%</td></tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
