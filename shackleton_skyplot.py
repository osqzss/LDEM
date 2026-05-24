"""
shackleton_skyplot.py

Draw a skyplot of the Shackleton crater-rim elevation mask from the CSV
produced by shackleton_edge_mask.c.

Skyplot convention:
  - Centre   : zenith (elevation = 90 deg)
  - Rim      : horizon (elevation = 0 deg)
  - Up       : North  (azimuth = 0 deg)
  - Right    : East   (azimuth = 90 deg)
  - Elevation grid   : every 30 deg
  - Azimuth grid     : every 30 deg (dotted lines)
  - Shaded region    : between the elevation mask and the horizon (blocked sky)

Usage:
  python shackleton_skyplot.py [input.csv] [output.png]

  Defaults:
    input  : shackleton_edge_mask.csv
    output : shackleton_skyplot.png
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import matplotlib.patheffects as pe

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

INPUT_CSV  = sys.argv[1] if len(sys.argv) > 1 else "shackleton_edge_mask.csv"
OUTPUT_PNG = sys.argv[2] if len(sys.argv) > 2 else "shackleton_skyplot.png"

ELEV_GRID_DEG = 30      # elevation grid interval [deg]
AZ_GRID_DEG   = 30      # azimuth  grid interval [deg]

FIG_SIZE   = (7, 7)
DPI        = 150

# Colours
COLOR_GRID        = "#888888"
COLOR_AZ_GRID     = "#aaaaaa"
COLOR_MASK_LINE   = "#c0392b"   # crater rim line
COLOR_FILL        = "#e74c3c"   # blocked-sky fill
FILL_ALPHA        = 0.20
COLOR_HORIZON     = "#555555"
COLOR_LABEL       = "#333333"

# ---------------------------------------------------------------------------
# Skyplot coordinate helper
# ---------------------------------------------------------------------------

def skyplot_xy(az_deg, el_deg):
    """
    Convert azimuth / elevation to Cartesian skyplot coordinates.

    Mapping:
      r = 1 - el_deg / 90   (zenith -> r=0, horizon -> r=1)
      theta = 0 at North, increases clockwise (East = 90 deg)

    Returns:
      x, y  (float or ndarray)
    """
    az  = np.deg2rad(az_deg)
    r   = 1.0 - el_deg / 90.0
    x   =  r * np.sin(az)   # East positive
    y   =  r * np.cos(az)   # North positive
    return x, y

# ---------------------------------------------------------------------------
# Load CSV
# ---------------------------------------------------------------------------

data = np.loadtxt(INPUT_CSV, delimiter=",", skiprows=1)
az_data  = data[:, 0]   # azimuth [deg]
el_data  = data[:, 1]   # elevation mask [deg]

# Close the polygon by appending the first point at the end.
az_closed = np.append(az_data, az_data[0])
el_closed = np.append(el_data, el_data[0])

mask_x, mask_y = skyplot_xy(az_closed, el_closed)

# Horizon circle for the fill polygon (closed, at el = 0).
az_circle = np.linspace(0.0, 360.0, 361)
hx, hy    = skyplot_xy(az_circle, np.zeros_like(az_circle))

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
ax.set_aspect("equal")
ax.axis("off")

# ---------------------------------------------------------------------------
# Elevation grid circles (every 30 deg: 0, 30, 60, 90)
# ---------------------------------------------------------------------------

for el in range(0, 91, ELEV_GRID_DEG):
    r = 1.0 - el / 90.0
    if r < 1e-9:
        # Zenith dot
        ax.plot(0, 0, "o", color=COLOR_GRID, markersize=2, zorder=3)
        continue
    theta = np.linspace(0, 2 * np.pi, 361)
    cx, cy = r * np.cos(theta), r * np.sin(theta)
    ax.plot(cx, cy,
            color=COLOR_GRID if el > 0 else COLOR_HORIZON,
            linewidth=1.2 if el == 0 else 0.7,
            linestyle="solid" if el == 0 else (0, (3, 4)),
            zorder=3)

    # Elevation label (placed at the North spoke, slightly right of it)
    if el > 0:
        lx = r * np.sin(np.deg2rad(2))  # tiny eastward offset
        ly = r * np.cos(np.deg2rad(2))
        ax.text(lx + 0.02, ly, f"{el}°",
                fontsize=7.5, color=COLOR_LABEL,
                va="center", ha="left", zorder=5)

# ---------------------------------------------------------------------------
# Azimuth grid lines (every 30 deg, dotted)
# ---------------------------------------------------------------------------

for az in range(0, 360, AZ_GRID_DEG):
    ex = np.sin(np.deg2rad(az))
    ey = np.cos(np.deg2rad(az))
    ax.plot([0, ex], [0, ey],
            color=COLOR_AZ_GRID, linewidth=0.7,
            linestyle=(0, (3, 4)),   # dotted
            zorder=2)

# ---------------------------------------------------------------------------
# Compass labels (N / E / S / W)
# ---------------------------------------------------------------------------

compass = {"N": (0, 1), "E": (1, 0), "S": (0, -1), "W": (-1, 0)}
for label, (lx, ly) in compass.items():
    ax.text(lx * 1.08, ly * 1.08, label,
            fontsize=11, fontweight="bold",
            color=COLOR_LABEL,
            ha="center", va="center", zorder=5)

# ---------------------------------------------------------------------------
# Blocked-sky fill (between elevation mask and horizon)
#
# Build a compound path:
#   outer contour = horizon circle (r=1, counter-clockwise)
#   inner contour = elevation mask polygon (clockwise)
# The even-odd fill rule then shades only the annular region.
# ---------------------------------------------------------------------------

# Horizon: counter-clockwise (azimuth 0 -> 360)
n_h    = len(az_circle)
codes_h = ([Path.MOVETO] + [Path.LINETO] * (n_h - 2) + [Path.CLOSEPOLY])
verts_h = list(zip(hx, hy))

# Elevation mask: clockwise (azimuth 360 -> 0) so that the winding rules
# create a hole inside the outer horizon polygon.
n_m    = len(mask_x)
codes_m = ([Path.MOVETO] + [Path.LINETO] * (n_m - 2) + [Path.CLOSEPOLY])
verts_m = list(zip(mask_x[::-1], mask_y[::-1]))

compound_path = Path(verts_h + verts_m,
                     codes_h + codes_m)

fill_patch = mpatches.PathPatch(compound_path,
                                facecolor=COLOR_FILL,
                                alpha=FILL_ALPHA,
                                edgecolor="none",
                                zorder=4)
ax.add_patch(fill_patch)

# ---------------------------------------------------------------------------
# Elevation mask line
# ---------------------------------------------------------------------------

ax.plot(mask_x, mask_y,
        color=COLOR_MASK_LINE,
        linewidth=1.8,
        zorder=5,
        label="Crater-rim elevation mask")

# ---------------------------------------------------------------------------
# Axes limits and title
# ---------------------------------------------------------------------------

margin = 1.18
ax.set_xlim(-margin, margin)
ax.set_ylim(-margin, margin)

ax.set_title("Shackleton Crater — Elevation Mask\n"
             r"$\phi_{rx}$ = −89.67°,  $\lambda_{rx}$ = 129.78°,  "
             r"$h_{ant}$ = 2.0 m",
             fontsize=11, pad=10)

ax.legend(loc="lower right", fontsize=8, framealpha=0.85)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches="tight")
print(f"Saved: {OUTPUT_PNG}")
