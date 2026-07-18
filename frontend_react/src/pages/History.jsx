import { useEffect, useState } from 'react';
import { getRun, listRuns } from '../api/client.js';
import WaferCanvas from '../components/WaferCanvas.jsx';

export default function History() {
  const [runs, setRuns] = useState(null);
  const [error, setError] = useState(null);
  const [selected, setSelected] = useState(null);

  useEffect(() => {
    listRuns()
      .then(setRuns)
      .catch((e) => setError(e.message));
  }, []);

  async function loadDetail(id) {
    try {
      const detail = await getRun(id);
      setSelected(detail);
    } catch (e) {
      setError(e.message);
    }
  }

  return (
    <div className="page">
      <h1>Run history</h1>
      {error && <p className="error-text">{error}</p>}
      {!runs && !error && <p className="muted">Loading\u2026</p>}
      {runs && runs.length === 0 && <p className="muted">No runs yet \u2014 go to Simulate and run a recipe.</p>}

      <div className="history-grid">
        <ul className="run-list">
          {runs?.map((r) => (
            <li key={r.id}>
              <button onClick={() => loadDetail(r.id)}>
                Run #{r.id} \u2014 {r.nx}x{r.ny} \u2014 {new Date(r.created_at).toLocaleString()}
              </button>
            </li>
          ))}
        </ul>
        {selected && (
          <div className="panel">
            <h3>Run #{selected.id}</h3>
            <WaferCanvas wafer={selected.wafer} />
          </div>
        )}
      </div>
    </div>
  );
}
