import { NavLink, Outlet } from 'react-router-dom';
import Logo from './Logo.jsx';

const links = [
  { to: '/', label: 'Simulate', end: true },
  { to: '/history', label: 'History' },
  { to: '/analytics', label: 'Analytics' },
];

export default function Layout() {
  return (
    <div className="app-shell">
      <aside className="sidebar">
        <div className="brand">
          <Logo size={28} />
          <span>AttoFab</span>
        </div>
        <nav>
          {links.map((l) => (
            <NavLink key={l.to} to={l.to} end={l.end} className={({ isActive }) => (isActive ? 'active' : '')}>
              {l.label}
            </NavLink>
          ))}
        </nav>
        <div className="sidebar-footer">Open-source fabrication simulator</div>
      </aside>
      <main className="content">
        <Outlet />
      </main>
    </div>
  );
}
