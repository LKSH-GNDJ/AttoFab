import { HashRouter, Routes, Route } from 'react-router-dom';
import Layout from './components/Layout.jsx';
import Simulate from './pages/Simulate.jsx';
import History from './pages/History.jsx';
import Analytics from './pages/Analytics.jsx';

export default function App() {
  return (
    <HashRouter>
      <Routes>
        <Route element={<Layout />}>
          <Route path="/" element={<Simulate />} />
          <Route path="/history" element={<History />} />
          <Route path="/analytics" element={<Analytics />} />
        </Route>
      </Routes>
    </HashRouter>
  );
}
