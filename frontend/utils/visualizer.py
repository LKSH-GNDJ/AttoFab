"""
utils/visualizer.py - renders a Wafer2d JSON export as a matplotlib
cross-section image (the Streamlit-native counterpart to web/index.html's
canvas renderer).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

MATERIAL_COLORS = {
    "Silicon": "#4a4e57",
    "Oxide": "#7fc8f8",
    "Photoresist": "#e8b84b",
    "PhotoresistExposed": "#f2d98a",
    "Metal": "#d4d8dc",
    "Void": "#0e1116",
}
MATERIAL_ORDER = list(MATERIAL_COLORS.keys())


def wafer_to_rgb_array(wafer: dict) -> np.ndarray:
    nx, ny = wafer["nx"], wafer["ny"]
    material = wafer["material"]
    idx = {m: i for i, m in enumerate(MATERIAL_ORDER)}
    palette = np.array([_hex_to_rgb(MATERIAL_COLORS[m]) for m in MATERIAL_ORDER])

    indices = np.array([idx[m] for m in material], dtype=np.int32).reshape(ny, nx)
    return palette[indices]


def _hex_to_rgb(hex_color: str) -> tuple[float, float, float]:
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def render_cross_section(wafer: dict, title: str | None = None):
    """Returns a matplotlib Figure showing the wafer cross-section,
    suitable for st.pyplot()."""
    rgb = wafer_to_rgb_array(wafer)
    nx, ny = wafer["nx"], wafer["ny"]
    dx_um, dy_um = wafer["dx_um"], wafer["dy_um"]

    fig, ax = plt.subplots(figsize=(8, 8 * ny / max(nx, 1)))
    ax.imshow(rgb, extent=[0, nx * dx_um, ny * dy_um, 0], aspect="auto", interpolation="nearest")
    ax.set_xlabel("x (\u00b5m)")
    ax.set_ylabel("depth (\u00b5m)")
    if title:
        ax.set_title(title)

    # Legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=MATERIAL_COLORS[m]) for m in MATERIAL_ORDER]
    ax.legend(handles, MATERIAL_ORDER, loc="upper right", fontsize=8, framealpha=0.85)

    fig.tight_layout()
    return fig
