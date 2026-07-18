import { useEffect, useRef } from 'react';

const MATERIAL_COLORS = {
  Silicon: '#4a4e57',
  Oxide: '#7fc8f8',
  Photoresist: '#e8b84b',
  PhotoresistExposed: '#f2d98a',
  Metal: '#d4d8dc',
  Void: '#0e1116',
};

export default function WaferCanvas({ wafer }) {
  const canvasRef = useRef(null);

  useEffect(() => {
    if (!wafer) return;
    const { nx, ny, material } = wafer;
    const canvas = canvasRef.current;
    const px = Math.max(2, Math.min(8, Math.floor(700 / nx)));
    canvas.width = nx * px;
    canvas.height = ny * px;

    const ctx = canvas.getContext('2d');
    const offscreen = document.createElement('canvas');
    offscreen.width = nx;
    offscreen.height = ny;
    const octx = offscreen.getContext('2d');
    const img = octx.createImageData(nx, ny);

    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const m = material[y * nx + x];
        const hex = MATERIAL_COLORS[m] || '#ff00ff';
        const r = parseInt(hex.slice(1, 3), 16);
        const g = parseInt(hex.slice(3, 5), 16);
        const b = parseInt(hex.slice(5, 7), 16);
        const i = (y * nx + x) * 4;
        img.data[i] = r;
        img.data[i + 1] = g;
        img.data[i + 2] = b;
        img.data[i + 3] = 255;
      }
    }
    octx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(offscreen, 0, 0, canvas.width, canvas.height);
  }, [wafer]);

  if (!wafer) return null;

  return (
    <div className="wafer-canvas-wrap">
      <canvas ref={canvasRef} className="wafer-canvas" />
      <div className="legend">
        {Object.entries(MATERIAL_COLORS).map(([name, color]) => (
          <span key={name} className="legend-item">
            <span className="swatch" style={{ background: color }} />
            {name}
          </span>
        ))}
      </div>
    </div>
  );
}
