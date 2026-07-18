import { useState } from 'react';
import { simulate } from '../api/client.js';
import WaferCanvas from '../components/WaferCanvas.jsx';
import ResultCard from '../components/ResultCard.jsx';

const DEFAULT_RECIPE = {
  nx: 60,
  ny: 80,
  dx_um: 0.01,
  dy_um: 0.01,
  substrate: { dopant: 'Boron', concentration_cm3: 1e15 },
  steps: [
    { op: 'oxidize', temperature_c: 1000, time_hours: 0.75, ambient: 'Dry' },
    { op: 'implant', dopant: 'Phosphorus', dose_cm2: 1e15, energy_kev: 80 },
    { op: 'anneal', temperature_c: 1000, time_minutes: 20 },
  ],
};

export default function Simulate() {
  const [recipeText, setRecipeText] = useState(JSON.stringify(DEFAULT_RECIPE, null, 2));
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null);

  async function run() {
    setError(null);
    setResult(null);
    let recipe;
    try {
      recipe = JSON.parse(recipeText);
    } catch (e) {
      setError(`Invalid JSON: ${e.message}`);
      return;
    }
    setLoading(true);
    try {
      const res = await simulate(recipe);
      setResult(res);
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  }

  return (
    <div className="page simulate-page">
      <h1>Simulate</h1>
      <div className="simulate-grid">
        <div className="panel">
          <h3>Process recipe</h3>
          <textarea
            value={recipeText}
            onChange={(e) => setRecipeText(e.target.value)}
            spellCheck={false}
            rows={24}
          />
          <button onClick={run} disabled={loading}>
            {loading ? 'Running\u2026' : 'Run simulation'}
          </button>
          {error && <p className="error-text">{error}</p>}
        </div>

        <div className="panel">
          <h3>Cross-section</h3>
          {result ? <WaferCanvas wafer={result.wafer} /> : <p className="muted">Run a recipe to see the wafer cross-section.</p>}
        </div>

        <div className="panel">
          <h3>Summary</h3>
          {result ? <ResultCard result={result} /> : <p className="muted">No result yet.</p>}
        </div>
      </div>
    </div>
  );
}
