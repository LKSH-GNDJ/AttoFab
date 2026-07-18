const BASE_URL = import.meta.env.VITE_ATTOFAB_API_BASE_URL || '/api';

async function request(path, options = {}) {
  const res = await fetch(`${BASE_URL}${path}`, {
    headers: { 'Content-Type': 'application/json' },
    ...options,
  });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`${res.status} ${res.statusText}: ${body}`);
  }
  return res.json();
}

export const health = () => request('/health');

export const simulate = (recipe) =>
  request('/simulate', { method: 'POST', body: JSON.stringify(recipe) });

export const listRuns = (limit = 50) => request(`/runs?limit=${limit}`);

export const getRun = (id) => request(`/runs/${id}`);
