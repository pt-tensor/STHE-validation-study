#!/usr/bin/env python3
"""Lab-scale SG-STHX structured blockMesh (Zahid et al., Thermal Science 2023).

Geometry from Table 1: 9 tubes (3×3), shell ID 44 mm, L_eff 186 mm, 4 segmental
baffles @ 22% cut / 35.6 mm spacing. Cylinders are O-grids.

Regions:
  hot_fluid  — tube ID, headers, axial ports
  cold_fluid — shell gaps, baffle holes/windows, nozzle stubs
  solid      — tube walls, baffles, sheets, shell wall, headers, nozzle walls

Run on the host:
  python3 scripts/makeStructuredBlockMesh.py
"""
from __future__ import print_function

import math
import os

# ---------------------------------------------------------------------------
# Geometry (m) — Zahid et al. Thermal Science 2023 Table 1 (SG-STHX lab-scale)
# DOI: 10.2298/TSCI220124075Z
# ---------------------------------------------------------------------------
# Tube OD = 6 mm, ID = 4 mm → wall = 1 mm
D_TUBE = 0.006
TUBE_WALL = 0.001
R_IN = 0.5 * D_TUBE - TUBE_WALL
R_OUT = 0.5 * D_TUBE
# Shell ID = 44 mm; OD = 50 mm → wall = 3 mm (paper: OD 50 mm, thickness 3 mm)
D_SHELL = 0.044
R_SHELL = 0.5 * D_SHELL
SHELL_WALL = 0.003
R_SHELL_OUT = R_SHELL + SHELL_WALL
# 3×3 square bank (9 tubes). Pitch not in Table 1; TEMA-like ~d+6 mm → 12 mm
S_L = 0.012
S_T = 0.012
N_SIDE = 3
N_TUBES = 9
L = 0.186                 # tube effective length (Table 1)
B_F = 0.0356              # baffle spacing
N_BAFFLES = 4
BAFFLE_THICK = 0.001
BAFFLE_CUT = 0.22
# Zero diametral tube–baffle hole clearance (tubes fit the baffle)
BAFFLE_HOLE_CLEARANCE = 0.0
MESH_HOLE_CLEARANCE = 0.0
R_HOLE = 0.5 * D_TUBE + 0.5 * MESH_HOLE_CLEARANCE  # == R_OUT
HOLE_COLLAPSED = abs(R_HOLE - R_OUT) < 1e-12

# Paper total length 200 mm = L_eff + 2×tube-sheet → sheet = 7 mm each
TUBE_SHEET = 0.007
HEADER_LEN = 0.010          # short plenums for tube ports (not in Table 1)
HEADER_WALL = 0.003
HEADER_MARGIN = 0.002
PORT_LEN = 0.010
# Fig. measure: D_shell / D_noz_eq = 8 → Ø5.5 mm circle.
# Rectangular opening: Z-length fixed at D_NOZZLE_EQ; X-width so W×L = A_circle.
D_NOZZLE_EQ = D_SHELL / 8.0
A_NOZZLE = math.pi * (0.5 * D_NOZZLE_EQ) ** 2
L_NOZZLE = D_NOZZLE_EQ              # Z opening length = 5.5 mm
W_NOZZLE = A_NOZZLE / L_NOZZLE      # X opening width ≈ 4.320 mm
R_NOZZLE = 0.5 * L_NOZZLE           # half Z-extent (nozzle_active / z-planes)
HALF_X_NOZZLE = 0.5 * W_NOZZLE
D_NOZZLE = W_NOZZLE
SIDE_NOZZLE = L_NOZZLE
NOZZLE_LEN = 0.020
NOZZLE_WALL = SHELL_WALL
TUBE_LAYOUT = "inline"      # 3×3 square
# Minimum radial thickness targets (scaled to lab shell)
MIN_PACK_ANNULUS = 0.002
MIN_SHELL_WALL = SHELL_WALL

# Mesh density — last good validation (~762k; Q≈1886 W / Δp≈1673 Pa).
# Taken from working-tree state immediately before b9b44544 (+25% / 1.52M).
# Selective Z from that commit; densities are the prior 762k baseline.
# N_CORE must equal N_THETA for conformal O-grids
N_CORE = 8
N_RAD_HOT = 3
N_WALL = 3          # ≥3 cells through tube wall thickness (1 mm)
N_RAD_HOLE = 1
N_RAD_ISO = 4          # uniform OD → iso circle (keeps OD corners unskewed)
N_RAD_COLD = 6         # iso → unit-cell morph
N_RAD_SHELL = 4          # pack → shell (cold) denser
N_RAD_CORNER = 3         # corner-tube exterior collar mid-split
N_SHELL_WALL = 1
N_THETA = 8
# Axial (Z): refine large shell gaps only; keep baffle + cold-nozzle slabs.
N_Z_GAP = 24              # large shell compartments (between baffles / away from nozzles)
N_Z_NOZZLE = 12           # keep density in cold inlet/outlet Z footprint
N_Z_BAFFLE = 3            # baffle slabs (slightly thicker resolution)
N_Z_SHEET = 4
N_Z_HEADER = 6
N_Z_PORT = 4
N_NOZZLE_Y = 8
BLEND_TS = (0.35, 0.7)    # two cold/hot blend rings (keep 2 for conformal stability)
# simpleGrading ratios (last/first cell size along block local axes)
G_RAD_ISO = 1.0           # uniform annulus at tube OD
G_RAD_COLD = 3.0          # iso → cell (milder; OD grading handled by iso layer)
G_RAD_HOT = 0.33          # core → rin: fine at tube ID
G_RAD_SHELL = 0.25        # pack → shell ID: fine at shell
G_RAD_CORNER = 1.0        # exterior collar mid-split (uniform)
G_SHELL_WALL = 1.0
G_Z_BAFFLE = 1.0
# Near-OD iso annulus (same on every tube). Sized from the tightest
# corner-tube clearance so mid/centre tubes match corner resolution.
MIN_CORNER_SHELL_GAP = 0.001
ISO_ANNULUS_FRAC = 0.55   # fraction of (r_iso_max − R_OUT) on corner tubes



# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
ANGLES = (
    -0.75 * math.pi,
    -0.25 * math.pi,
    0.25 * math.pi,
    0.75 * math.pi,
)


def _is_pack_perimeter(row, col, k):
    """True if unit-cell corner k lies on the outer bank perimeter.

    k: 0=SW, 1=SE, 2=NE, 3=NW of that tube's square.
    """
    if row == 0 and k in (0, 1):
        return True
    if row == N_SIDE - 1 and k in (2, 3):
        return True
    if col == 0 and k in (0, 3):
        return True
    if col == N_SIDE - 1 and k in (1, 2):
        return True
    return False


def _outer_diagonal_k(row, col):
    """Outer diagonal corner index for a corner tube, else None.

    That unique corner is not shared with another tube and is the site of
    the tight tube–shell gap along the bank diagonal.
    """
    if row == 0 and col == 0:
        return 0  # SW
    if row == 0 and col == N_SIDE - 1:
        return 1  # SE
    if row == N_SIDE - 1 and col == N_SIDE - 1:
        return 2  # NE
    if row == N_SIDE - 1 and col == 0:
        return 3  # NW
    return None


def _circularize_xy(x, y, radius):
    r = math.hypot(x, y)
    if r < 1e-14:
        return (radius, 0.0)
    return (radius * x / r, radius * y / r)


def _lerp_xy(a, b, t):
    return (a[0] * (1.0 - t) + b[0] * t, a[1] * (1.0 - t) + b[1] * t)


def _blend_square_circle(cx, cy, half, radius, t):
    """Corner ring morphing square(half) → circle(radius). t=0 square, t=1 circle."""
    sq = _square_xy(cx, cy, half)
    circ = _ring_xy(cx, cy, radius)
    return [_lerp_xy(sq[i], circ[i], t) for i in range(4)]


def _blend_arc_mid(cx, cy, half, radius, k, t):
    """Edge midpoint morphing straight square edge → circular arc."""
    sq = _square_xy(cx, cy, half)
    k2 = (k + 1) % 4
    sm = (0.5 * (sq[k][0] + sq[k2][0]), 0.5 * (sq[k][1] + sq[k2][1]))
    cm = _mid_arc(cx, cy, radius, k)
    return _lerp_xy(sm, cm, t)


def tube_center_xy():
    """Return N_TUBES (x, y) centres [m].

    inline  → square N_SIDE×N_SIDE lattice (Zahid lab: 3×3).
    staggered → legacy 4×4 STHE layout (kept for older CAD cases).
    """
    pts = []
    if TUBE_LAYOUT == "inline":
        # Centres at …, -S, 0, +S, … with S = pitch
        half = 0.5 * (N_SIDE - 1)
        for j in range(N_SIDE):
            y = (j - half) * S_T
            for i in range(N_SIDE):
                x = (i - half) * S_L
                pts.append((x, y))
    else:
        y_rows = [-1.5 * S_T, -0.5 * S_T, 0.5 * S_T, 1.5 * S_T]
        for i, y in enumerate(y_rows):
            if i % 2 == 0:
                xs = [-1.5 * S_L, -0.5 * S_L, 0.5 * S_L, 1.5 * S_L]
            else:
                xs = [-1.0 * S_L, 0.0, 1.0 * S_L, 2.0 * S_L]
                xs = [x - 0.5 * S_L for x in xs]
            for x in xs:
                pts.append((x, y))
    if len(pts) != N_TUBES:
        raise RuntimeError("Expected %d tubes, got %d" % (N_TUBES, len(pts)))
    r_max = R_SHELL - R_OUT
    for x, y in pts:
        if math.hypot(x, y) > r_max + 1e-9:
            raise RuntimeError("Tube (%.4f,%.4f) outside shell clearance" % (x, y))
    return pts


def header_inner_radius(centers):
    pack = max(math.hypot(x, y) + R_OUT for x, y in centers)
    return max(R_SHELL, pack + HEADER_MARGIN)


def _jac(p0, p1, p3, p4):
    e1 = (p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2])
    e2 = (p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2])
    e3 = (p4[0] - p0[0], p4[1] - p0[1], p4[2] - p0[2])
    return (
        e1[0] * (e2[1] * e3[2] - e2[2] * e3[1])
        + e1[1] * (e2[2] * e3[0] - e2[0] * e3[2])
        + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0])
    )


def _ring_xy(cx, cy, radius):
    return [
        (cx + radius * math.cos(a), cy + radius * math.sin(a))
        for a in ANGLES
    ]


def _square_xy(cx, cy, half):
    return [
        (cx - half, cy - half),
        (cx + half, cy - half),
        (cx + half, cy + half),
        (cx - half, cy + half),
    ]


def _mid_arc(cx, cy, radius, k):
    k2 = (k + 1) % 4
    span = (ANGLES[k2] - ANGLES[k]) % (2.0 * math.pi)
    mid_a = ANGLES[k] + 0.5 * span
    return (cx + radius * math.cos(mid_a), cy + radius * math.sin(mid_a))


def _ring_xz(cx, cz, radius):
    """Circle in XZ plane (for nozzles along ±Y). Order matches ANGLES about +Y."""
    # ANGLES in (x,z): cos→x, sin→z
    return [
        (cx + radius * math.cos(a), cz + radius * math.sin(a))
        for a in ANGLES
    ]


def _square_xz(cx, cz, half):
    return [
        (cx - half, cz - half),
        (cx + half, cz - half),
        (cx + half, cz + half),
        (cx - half, cz + half),
    ]


def _mid_arc_xz(cx, cz, radius, k):
    k2 = (k + 1) % 4
    span = (ANGLES[k2] - ANGLES[k]) % (2.0 * math.pi)
    mid_a = ANGLES[k] + 0.5 * span
    return (cx + radius * math.cos(mid_a), cz + radius * math.sin(mid_a))


class Mesh(object):
    def __init__(self):
        self.pts = []
        self._key = {}
        self.blocks = []
        self.arcs = {}
        self.faces = {
            "inlet_hot": [],
            "outlet_hot": [],
            "inlet_cold": [],
            "outlet_cold": [],
            "shellOuter": [],
            "defaultWall": [],
        }

    def add(self, x, y, z):
        key = (round(x, 10), round(y, 10), round(z, 10))
        if key in self._key:
            return self._key[key]
        i = len(self.pts)
        self.pts.append((float(x), float(y), float(z)))
        self._key[key] = i
        return i

    def arc(self, i, j, mid):
        a, b = (i, j) if i < j else (j, i)
        self.arcs[(a, b)] = mid

    def block(self, zone, v, div, face_map=None, grading=(1.0, 1.0, 1.0)):
        # Keep (nx, ny) aligned with edges 0-1 and 0-3. If winding is flipped
        # to fix Jacobian, those edges swap — so swap nx/ny too. Otherwise
        # adjacent blocks disagree on face subdivisions whenever nx != ny.
        grading = (float(grading[0]), float(grading[1]), float(grading[2]))
        if _jac(self.pts[v[0]], self.pts[v[1]], self.pts[v[3]], self.pts[v[4]]) <= 0.0:
            v = (v[0], v[3], v[2], v[1], v[4], v[7], v[6], v[5])
            div = (div[1], div[0], div[2])
            grading = (grading[1], grading[0], grading[2])
            if face_map:
                # Face vertex order follows the swapped hex corners.
                remapped = []
                for name, face in face_map:
                    remapped.append((name, (face[0], face[3], face[2], face[1])))
                face_map = remapped
            if _jac(self.pts[v[0]], self.pts[v[1]], self.pts[v[3]], self.pts[v[4]]) <= 0.0:
                raise RuntimeError("inverted block in zone %s verts=%s" % (zone, v))
        self.blocks.append((zone, v, div, grading))
        if face_map:
            for name, face in face_map:
                self.faces[name].append(face)


def build_z_planes():
    zs = set()
    z_port_in = -(TUBE_SHEET + HEADER_LEN + HEADER_WALL + PORT_LEN)
    z_head_in_wall0 = -(TUBE_SHEET + HEADER_LEN + HEADER_WALL)
    z_head_in0 = -(TUBE_SHEET + HEADER_LEN)
    z_sheet0 = -TUBE_SHEET
    z0 = 0.0
    zL = L
    z_sheet1 = L + TUBE_SHEET
    z_head_out1 = L + TUBE_SHEET + HEADER_LEN
    z_head_out_wall1 = L + TUBE_SHEET + HEADER_LEN + HEADER_WALL
    z_port_out = L + TUBE_SHEET + HEADER_LEN + HEADER_WALL + PORT_LEN

    for z in (
        z_port_in, z_head_in_wall0, z_head_in0, z_sheet0, z0,
        zL, z_sheet1, z_head_out1, z_head_out_wall1, z_port_out,
    ):
        zs.add(round(z, 12))

    for i in range(1, N_BAFFLES + 1):
        zc = i * B_F
        zs.add(round(zc - 0.5 * BAFFLE_THICK, 12))
        zs.add(round(zc + 0.5 * BAFFLE_THICK, 12))

    for zc in (0.5 * B_F, L - 0.5 * B_F):
        zs.add(round(zc - R_NOZZLE, 12))
        zs.add(round(zc + R_NOZZLE, 12))

    return sorted(zs)


def layer_kind(z0, z1):
    zm = 0.5 * (z0 + z1)
    if zm < -(TUBE_SHEET + HEADER_LEN + HEADER_WALL):
        return "hot_port_in"
    if zm < -(TUBE_SHEET + HEADER_LEN):
        return "header_wall_in"
    if zm < -TUBE_SHEET:
        return "header_cav_in"
    if zm < 0.0:
        return "tube_sheet_in"
    if zm > L + TUBE_SHEET + HEADER_LEN + HEADER_WALL:
        return "hot_port_out"
    if zm > L + TUBE_SHEET + HEADER_LEN:
        return "header_wall_out"
    if zm > L + TUBE_SHEET:
        return "header_cav_out"
    if zm > L:
        return "tube_sheet_out"
    for i in range(1, N_BAFFLES + 1):
        zc = i * B_F
        if abs(zm - zc) <= 0.5 * BAFFLE_THICK + 1e-9:
            return "baffle", i
    return "shell_gap"


def baffle_chord_y(baffle_index):
    """Segmental baffle chord y = ±(R_SHELL − cut·D_SHELL).

    Odd baffles: open +Y (chord at +y_chord). Even: open −Y (chord at −y_chord).
    """
    y_chord = R_SHELL - BAFFLE_CUT * D_SHELL
    if baffle_index % 2 == 1:
        return y_chord
    return -y_chord


def baffle_in_window(x, y, baffle_index):
    """True if (x,y) is on the open (cold) side of the segmental chord."""
    yc = baffle_chord_y(baffle_index)
    # Orientation from 5b285521 / good 762k era (pre–b9b44544 flip).
    if baffle_index % 2 == 1:
        return y > yc  # open +Y
    return y < yc     # open -Y


def edge_intersect_chord_y(p0, p1, yc, eps=1e-14):
    """Intersection of segment p0→p1 with horizontal line y=yc, or None.

    Returns (x, yc, z) using linear interpolation in x and z.
    """
    y0, y1 = p0[1], p1[1]
    if abs(y1 - y0) < eps:
        return None
    if (y0 - yc) * (y1 - yc) > 0.0:
        return None
    t = (yc - y0) / (y1 - y0)
    if t < -1e-9 or t > 1.0 + 1e-9:
        return None
    t = max(0.0, min(1.0, t))
    x = p0[0] + t * (p1[0] - p0[0])
    z = p0[2] + t * (p1[2] - p0[2])
    return (x, yc, z, t)


def sector_straddle(corners, yc, eps=1e-12):
    """True if 2D/3D corner list has points on both sides of y=yc."""
    has_pos = False
    has_neg = False
    for c in corners:
        dy = c[1] - yc
        if dy > eps:
            has_pos = True
        elif dy < -eps:
            has_neg = True
    return has_pos and has_neg


def _theta_cell_split(t, n_total):
    """Split n_total cells across parameter t in (0,1); both sides ≥ 1."""
    n_total = max(2, int(n_total))
    n_a = int(round(t * n_total))
    n_a = max(1, min(n_total - 1, n_a))
    return n_a, n_total - n_a


def nz_for(kind, zm=None):
    if kind == "baffle":
        return N_Z_BAFFLE
    if kind in ("tube_sheet_in", "tube_sheet_out"):
        return N_Z_SHEET
    if kind in ("header_cav_in", "header_cav_out", "header_wall_in", "header_wall_out"):
        return N_Z_HEADER
    if kind in ("hot_port_in", "hot_port_out"):
        return N_Z_PORT
    # shell_gap: keep cold-nozzle Z bands; refine large compartments
    if zm is not None and (
        nozzle_active(zm, "inlet") or nozzle_active(zm, "outlet")
    ):
        return N_Z_NOZZLE
    return N_Z_GAP


def nozzle_active(zm, side):
    """Whether mid-z is within the nozzle axial footprint.

    Counterflow (paper): cold inlet near hot outlet (+Z), cold outlet near hot
    inlet (−Z / z≈0). Hot ports: inlet at −Z, outlet at +Z.
    """
    zc = (L - 0.5 * B_F) if side == "inlet" else (0.5 * B_F)
    return abs(zm - zc) <= R_NOZZLE + 1e-9


def build():
    centers = tube_center_xy()
    r_header = header_inner_radius(centers)
    zs = build_z_planes()
    m = Mesh()
    half_cell = 0.5 * S_L
    s_core = 0.50 * R_IN
    bank_half = 0.5 * N_SIDE * S_L
    # Pack circle: keep a usable annulus to the shell (no crushed corner gap)
    r_tube_extent = max(math.hypot(cx, cy) + R_OUT for cx, cy in centers)
    r_pack = max(
        r_tube_extent + 0.003,
        0.5 * (bank_half + bank_half * math.sqrt(2.0)),
    )
    r_pack = min(r_pack, R_SHELL - MIN_PACK_ANNULUS)
    if r_pack <= r_tube_extent + 0.001:
        r_pack = 0.5 * (r_tube_extent + R_SHELL)

    # Uniform near-OD iso annulus on EVERY tube (= corner-tube thickness).
    # Corner clearance sets the target so mid tubes are not coarser at the OD.
    iso_annulus = None
    for cx, cy in centers:
        r_iso_max = R_SHELL - MIN_CORNER_SHELL_GAP - math.hypot(cx, cy)
        span = r_iso_max - R_OUT
        if span > 1e-9:
            ann = ISO_ANNULUS_FRAC * span
            iso_annulus = ann if iso_annulus is None else min(iso_annulus, ann)
    if iso_annulus is None:
        iso_annulus = ISO_ANNULUS_FRAC * (0.72 * half_cell - R_OUT)
    # Keep iso inside the unit-cell square for the iso→cell morph.
    iso_annulus = min(iso_annulus, max(0.72 * half_cell - R_OUT, 1e-6))

    # Per-z rings: core / hot-blend / rin / rout / rhole / iso / cold-blend / cell
    # Iso = tube-centred circle just outside OD (uniform annulus → no OD corner crush).
    tube_core, tube_rin, tube_rout, tube_rhole, tube_iso, tube_cell = (
        [], [], [], [], [], []
    )
    tube_hot_b = []   # list of blend rings core→rin (len = len(BLEND_TS))
    tube_cold_b = []  # list of blend rings iso→cell
    shell_in, shell_out, bank_sq, header_r = [], [], [], []

    for z in zs:
        tc, tri, tro, trh, tiso, tcell = [], [], [], [], [], []
        thb = [[] for _ in BLEND_TS]
        tcb = [[] for _ in BLEND_TS]
        for it, (cx, cy) in enumerate(centers):
            row, col = divmod(it, N_SIDE)
            diag_k = _outer_diagonal_k(row, col)
            core = _square_xy(cx, cy, s_core)
            rin = _ring_xy(cx, cy, R_IN)
            rout = _ring_xy(cx, cy, R_OUT)
            rhole = _ring_xy(cx, cy, R_HOLE)
            r_iso = R_OUT + iso_annulus
            r_iso_max = R_SHELL - MIN_CORNER_SHELL_GAP - math.hypot(cx, cy)
            r_cell_diag = r_iso_max if diag_k is not None else None
            if r_iso <= R_OUT + 1e-6:
                raise RuntimeError(
                    "iso ring fits inside OD at tube (%.4f,%.4f)" % (cx, cy)
                )
            if diag_k is not None and r_iso >= r_iso_max - 1e-9:
                raise RuntimeError(
                    "iso ring reaches shell at corner tube (%.4f,%.4f)" % (cx, cy)
                )
            iso = _ring_xy(cx, cy, r_iso)

            cell = _square_xy(cx, cy, half_cell)
            for k in range(4):
                if _is_pack_perimeter(row, col, k):
                    if diag_k is not None and k == diag_k:
                        # Tube-radial outer corner (not origin-pack snap).
                        cell[k] = _ring_xy(cx, cy, r_cell_diag)[k]
                    else:
                        cell[k] = _circularize_xy(cell[k][0], cell[k][1], r_pack)

            # Smooth square→circle (core→rin) and iso→cell (cold outer morph)
            hot_blends = [
                _blend_square_circle(cx, cy, s_core, R_IN, t) for t in BLEND_TS
            ]
            cold_blends = []
            for t in BLEND_TS:
                cold_blends.append([
                    _lerp_xy(iso[k], cell[k], t) for k in range(4)
                ])

            tc.append([m.add(p[0], p[1], z) for p in core])
            for ib, ring in enumerate(hot_blends):
                thb[ib].append([m.add(p[0], p[1], z) for p in ring])
            tri.append([m.add(p[0], p[1], z) for p in rin])
            tro.append([m.add(p[0], p[1], z) for p in rout])
            if HOLE_COLLAPSED:
                # Reuse rout verts — no degenerate hole ring
                trh.append(tro[-1])
            else:
                trh.append([m.add(p[0], p[1], z) for p in rhole])
            tiso.append([m.add(p[0], p[1], z) for p in iso])
            for ib, ring in enumerate(cold_blends):
                tcb[ib].append([m.add(p[0], p[1], z) for p in ring])
            tcell.append([m.add(p[0], p[1], z) for p in cell])

            for k in range(4):
                k2 = (k + 1) % 4
                for ib, t in enumerate(BLEND_TS):
                    mx, my = _blend_arc_mid(cx, cy, s_core, R_IN, k, t)
                    m.arc(thb[ib][-1][k], thb[ib][-1][k2], (mx, my, z))
                mx, my = _mid_arc(cx, cy, R_IN, k)
                m.arc(tri[-1][k], tri[-1][k2], (mx, my, z))
                mx, my = _mid_arc(cx, cy, R_OUT, k)
                m.arc(tro[-1][k], tro[-1][k2], (mx, my, z))
                if not HOLE_COLLAPSED:
                    mx, my = _mid_arc(cx, cy, R_HOLE, k)
                    m.arc(trh[-1][k], trh[-1][k2], (mx, my, z))
                mx, my = _mid_arc(cx, cy, r_iso, k)
                m.arc(tiso[-1][k], tiso[-1][k2], (mx, my, z))
                # Cold blends: lerp iso-arc mid → cell-edge mid
                for ib, t in enumerate(BLEND_TS):
                    hx, hy = _mid_arc(cx, cy, r_iso, k)
                    cxm = 0.5 * (cell[k][0] + cell[k2][0])
                    cym = 0.5 * (cell[k][1] + cell[k2][1])
                    mx, my = _lerp_xy((hx, hy), (cxm, cym), t)
                    m.arc(tcb[ib][-1][k], tcb[ib][-1][k2], (mx, my, z))
        tube_core.append(tc)
        tube_hot_b.append(thb)
        tube_rin.append(tri)
        tube_rout.append(tro)
        tube_rhole.append(trh)
        tube_iso.append(tiso)
        tube_cold_b.append(tcb)
        tube_cell.append(tcell)

        si = _ring_xy(0.0, 0.0, R_SHELL)
        so = _ring_xy(0.0, 0.0, R_SHELL_OUT)
        shell_in.append([m.add(p[0], p[1], z) for p in si])
        shell_out.append([m.add(p[0], p[1], z) for p in so])
        for k in range(4):
            k2 = (k + 1) % 4
            mx, my = _mid_arc(0.0, 0.0, R_SHELL, k)
            m.arc(shell_in[-1][k], shell_in[-1][k2], (mx, my, z))
            mx, my = _mid_arc(0.0, 0.0, R_SHELL_OUT, k)
            m.arc(shell_out[-1][k], shell_out[-1][k2], (mx, my, z))

        bsq = _square_xy(0.0, 0.0, bank_half)
        bank_sq.append([m.add(p[0], p[1], z) for p in bsq])

        hr = _ring_xy(0.0, 0.0, r_header)
        header_r.append([m.add(p[0], p[1], z) for p in hr])
        for k in range(4):
            k2 = (k + 1) % 4
            mx, my = _mid_arc(0.0, 0.0, r_header, k)
            m.arc(header_r[-1][k], header_r[-1][k2], (mx, my, z))

    def add_prism(zone, bot, top, div_xy, nz, face_map=None, grading=(1.0, 1.0, 1.0)):
        v = (bot[0], bot[1], bot[2], bot[3], top[0], top[1], top[2], top[3])
        m.block(zone, v, (div_xy[0], div_xy[1], nz), face_map, grading)

    def add_ogrid_sector(zone, b_in, b_out, t_in, t_out, n_rad, nz, face_map=None,
                         grading=(1.0, 1.0, 1.0)):
        for k in range(4):
            k2 = (k + 1) % 4
            v = (
                b_in[k], b_out[k], b_out[k2], b_in[k2],
                t_in[k], t_out[k], t_out[k2], t_in[k2],
            )
            fmap = None
            if face_map:
                fmap = []
                for name, which in face_map:
                    if which == "outer":
                        fmap.append((name, (v[1], v[2], v[6], v[5])))
                    elif which == "bot":
                        fmap.append((name, (v[0], v[3], v[2], v[1])))
                    elif which == "top":
                        fmap.append((name, (v[4], v[5], v[6], v[7])))
                    elif which == ("sector", k):
                        fmap.append((name, (v[1], v[2], v[6], v[5])))
            m.block(zone, v, (n_rad, N_THETA, nz), fmap, grading)

    def add_ogrid_sector_k(zone, b_in, b_out, t_in, t_out, n_rad, nz, k, face_map=None,
                           grading=(1.0, 1.0, 1.0)):
        """Single O-grid sector k only."""
        k2 = (k + 1) % 4
        v = (
            b_in[k], b_out[k], b_out[k2], b_in[k2],
            t_in[k], t_out[k], t_out[k2], t_in[k2],
        )
        m.block(zone, v, (n_rad, N_THETA, nz), face_map, grading)

    def _split_radial(n_total, n_layers):
        base, rem = divmod(max(n_total, n_layers), n_layers)
        return [max(1, base + (1 if i < rem else 0)) for i in range(n_layers)]

    def add_morph_chain(zone, b_rings, t_rings, n_total, nz, face_map=None,
                        grading=(1.0, 1.0, 1.0)):
        """Connect concentric rings with square↔circle morph, splitting n_total."""
        n_layers = len(b_rings) - 1
        ns = _split_radial(n_total, n_layers)
        for i, nr in enumerate(ns):
            fmap = face_map if (face_map and i == 0) else None
            add_ogrid_sector(
                zone, b_rings[i], b_rings[i + 1],
                t_rings[i], t_rings[i + 1], nr, nz, fmap, grading,
            )

    def add_morph_chain_k(zone, b_rings, t_rings, n_total, nz, k,
                          grading=(1.0, 1.0, 1.0)):
        n_layers = len(b_rings) - 1
        ns = _split_radial(n_total, n_layers)
        for i, nr in enumerate(ns):
            add_ogrid_sector_k(
                zone, b_rings[i], b_rings[i + 1],
                t_rings[i], t_rings[i + 1], nr, nz, k, None, grading,
            )

    def add_cold_ogrid(zone, b_rout, t_rout, b_iso, t_iso, b_cold_outer, t_cold_outer,
                       nz, face_map=None, grading_outer=(1.0, 1.0, 1.0)):
        """OD → iso (uniform) then iso → cell morph (blends + cell)."""
        add_ogrid_sector(
            zone, b_rout, b_iso, t_rout, t_iso, N_RAD_ISO, nz, face_map,
            (G_RAD_ISO, 1.0, grading_outer[2]),
        )
        add_morph_chain(
            zone, b_cold_outer, t_cold_outer, N_RAD_COLD, nz, None, grading_outer,
        )

    def add_cold_ogrid_k(zone, b_rout, t_rout, b_iso, t_iso, b_cold_outer, t_cold_outer,
                         nz, k, grading_outer=(1.0, 1.0, 1.0)):
        add_ogrid_sector_k(
            zone, b_rout, b_iso, t_rout, t_iso, N_RAD_ISO, nz, k, None,
            (G_RAD_ISO, 1.0, grading_outer[2]),
        )
        add_morph_chain_k(
            zone, b_cold_outer, t_cold_outer, N_RAD_COLD, nz, k, grading_outer,
        )

    def _shell_chord_xy(xa, ya, xb, yb, radius, yc):
        """Point on circle radius at y=yc lying on directed arc a→b, else None."""
        disc = radius * radius - yc * yc
        if disc < 0.0:
            return None
        xs = math.sqrt(disc)
        a0 = math.atan2(ya, xa)
        a1 = math.atan2(yb, xb)
        da = (a1 - a0 + math.pi) % (2.0 * math.pi) - math.pi

        def _on_span(x, y):
            t = math.atan2(y, x)
            dt = (t - a0 + math.pi) % (2.0 * math.pi) - math.pi
            if da >= 0.0:
                return -1e-9 <= dt <= da + 1e-9
            return da - 1e-9 <= dt <= 1e-9

        for x in (xs, -xs):
            if _on_span(x, yc):
                return (x, yc)
        return None

    def add_ogrid_sector_k_baffle(
        b_in, b_out, t_in, t_out, n_rad, nz, k, grading=(1.0, 1.0, 1.0),
        outer_radius=None, force_zone=None, face_map_outer=None,
        theta_split=None,
    ):
        """O-grid sector on a baffle slab: full-zone or chord-split into 2 hexes.

        force_zone: if set, both pieces use this zone (e.g. shell OD metal).
        face_map_outer: optional patch name for the outer radial face (shellOuter).
        theta_split: optional (n_a, n_b) forcing θ cell counts for conformal neighbours.
        """
        k2 = (k + 1) % 4
        yc = baffle_chord_y(baff_i)
        p = m.pts
        bot = (b_in[k], b_out[k], b_out[k2], b_in[k2])
        top = (t_in[k], t_out[k], t_out[k2], t_in[k2])
        corners = [p[i] for i in bot]

        def _zone_for(v_h):
            if force_zone is not None:
                return force_zone
            mx = 0.25 * sum(m.pts[i][0] for i in v_h[:4])
            my = 0.25 * sum(m.pts[i][1] for i in v_h[:4])
            return "cold_fluid" if baffle_in_window(mx, my, baff_i) else "solid"

        def _fmap(v_h, zone):
            if not face_map_outer:
                return None
            return [(face_map_outer, (v_h[1], v_h[2], v_h[6], v_h[5]))]

        if not sector_straddle(corners, yc):
            mx = 0.25 * sum(p[i][0] for i in bot)
            my = 0.25 * sum(p[i][1] for i in bot)
            zone = force_zone if force_zone is not None else (
                "cold_fluid" if baffle_in_window(mx, my, baff_i) else "solid"
            )
            v = (
                b_in[k], b_out[k], b_out[k2], b_in[k2],
                t_in[k], t_out[k], t_out[k2], t_in[k2],
            )
            m.block(zone, v, (n_rad, N_THETA, nz), _fmap(v, zone), grading)
            return

        # Prefer θ-split (chord crosses inner 0–3 and outer 1–2 edges)
        hit_in_b = edge_intersect_chord_y(p[bot[0]], p[bot[3]], yc)
        hit_out_b = edge_intersect_chord_y(p[bot[1]], p[bot[2]], yc)
        hit_in_t = edge_intersect_chord_y(p[top[0]], p[top[3]], yc)
        hit_out_t = edge_intersect_chord_y(p[top[1]], p[top[2]], yc)

        if outer_radius is not None and hit_out_b is not None:
            xy = _shell_chord_xy(
                p[bot[1]][0], p[bot[1]][1], p[bot[2]][0], p[bot[2]][1],
                outer_radius, yc,
            )
            if xy is not None:
                hit_out_b = (xy[0], xy[1], p[bot[1]][2], hit_out_b[3])
                t_param = hit_out_t[3] if hit_out_t is not None else hit_out_b[3]
                hit_out_t = (xy[0], xy[1], p[top[1]][2], t_param)

        if (
            hit_in_b is not None and hit_out_b is not None
            and hit_in_t is not None and hit_out_t is not None
        ):
            f_b = m.add(hit_in_b[0], hit_in_b[1], hit_in_b[2])
            c_b = m.add(hit_out_b[0], hit_out_b[1], hit_out_b[2])
            f_t = m.add(hit_in_t[0], hit_in_t[1], hit_in_t[2])
            c_t = m.add(hit_out_t[0], hit_out_t[1], hit_out_t[2])
            if theta_split is not None:
                n_a, n_b = theta_split
            else:
                # Prefer outer-edge t so pack/cell shared edges match across blocks
                n_a, n_b = _theta_cell_split(hit_out_b[3], N_THETA)
            v_a = (
                bot[0], bot[1], c_b, f_b,
                top[0], top[1], c_t, f_t,
            )
            v_b = (
                f_b, c_b, bot[2], bot[3],
                f_t, c_t, top[2], top[3],
            )
            if outer_radius is not None:
                def _arc_pair(i0, i1, radius):
                    x0, y0, z = m.pts[i0][0], m.pts[i0][1], m.pts[i0][2]
                    x1, y1 = m.pts[i1][0], m.pts[i1][1]
                    a0 = math.atan2(y0, x0)
                    a1 = math.atan2(y1, x1)
                    da = (a1 - a0 + math.pi) % (2.0 * math.pi) - math.pi
                    am = a0 + 0.5 * da
                    m.arc(i0, i1, (radius * math.cos(am), radius * math.sin(am), z))
                _arc_pair(bot[1], c_b, outer_radius)
                _arc_pair(c_b, bot[2], outer_radius)
                _arc_pair(top[1], c_t, outer_radius)
                _arc_pair(c_t, top[2], outer_radius)
            for v_h, n_th in ((v_a, n_a), (v_b, n_b)):
                zone = _zone_for(v_h)
                m.block(zone, v_h, (n_rad, n_th, nz), _fmap(v_h, zone), grading)
            return

        # Radial split: chord crosses edges 0–1 and 3–2
        hit_r0_b = edge_intersect_chord_y(p[bot[0]], p[bot[1]], yc)
        hit_r1_b = edge_intersect_chord_y(p[bot[3]], p[bot[2]], yc)
        hit_r0_t = edge_intersect_chord_y(p[top[0]], p[top[1]], yc)
        hit_r1_t = edge_intersect_chord_y(p[top[3]], p[top[2]], yc)
        if (
            hit_r0_b is not None and hit_r1_b is not None
            and hit_r0_t is not None and hit_r1_t is not None
        ):
            c0_b = m.add(hit_r0_b[0], hit_r0_b[1], hit_r0_b[2])
            c1_b = m.add(hit_r1_b[0], hit_r1_b[1], hit_r1_b[2])
            c0_t = m.add(hit_r0_t[0], hit_r0_t[1], hit_r0_t[2])
            c1_t = m.add(hit_r1_t[0], hit_r1_t[1], hit_r1_t[2])
            n_a, n_b = _theta_cell_split(0.5 * (hit_r0_b[3] + hit_r1_b[3]), n_rad)
            v_in = (
                bot[0], c0_b, c1_b, bot[3],
                top[0], c0_t, c1_t, top[3],
            )
            v_out = (
                c0_b, bot[1], bot[2], c1_b,
                c0_t, top[1], top[2], c1_t,
            )
            for v_h, nr, g in (
                (v_in, n_a, grading),
                (v_out, n_b, (1.0, grading[1], grading[2])),
            ):
                zone = _zone_for(v_h)
                m.block(zone, v_h, (nr, N_THETA, nz), _fmap(v_h, zone), g)
            return

        # Fallback: centroid zone
        mx = 0.25 * sum(p[i][0] for i in bot)
        my = 0.25 * sum(p[i][1] for i in bot)
        zone = force_zone if force_zone is not None else (
            "cold_fluid" if baffle_in_window(mx, my, baff_i) else "solid"
        )
        v = (
            b_in[k], b_out[k], b_out[k2], b_in[k2],
            t_in[k], t_out[k], t_out[k2], t_in[k2],
        )
        m.block(zone, v, (n_rad, N_THETA, nz), _fmap(v, zone), grading)

    def add_morph_chain_k_baffle(b_rings, t_rings, n_total, nz, k,
                                grading=(1.0, 1.0, 1.0)):
        n_layers = len(b_rings) - 1
        ns = _split_radial(n_total, n_layers)
        # Canonical θ-split from outermost ring edge (matches pack→shell collar)
        yc = baffle_chord_y(baff_i)
        k2 = (k + 1) % 4
        outer = b_rings[-1]
        hit = edge_intersect_chord_y(
            m.pts[outer[k]], m.pts[outer[k2]], yc
        )
        theta_split = _theta_cell_split(hit[3], N_THETA) if hit is not None else None
        for i, nr in enumerate(ns):
            add_ogrid_sector_k_baffle(
                b_rings[i], b_rings[i + 1],
                t_rings[i], t_rings[i + 1],
                nr, nz, k, grading, theta_split=theta_split,
            )

    def add_cold_ogrid_k_baffle(
        b_rout, t_rout, b_iso, t_iso, b_outer, t_outer, nz, k,
        grading_outer=(1.0, 1.0, 1.0),
    ):
        # One θ-split from the outermost cell edge so every radial layer
        # (OD→iso and iso→cell) agrees with the pack→shell collar.
        yc = baffle_chord_y(baff_i)
        k2 = (k + 1) % 4
        outer = b_outer[-1]
        hit = edge_intersect_chord_y(
            m.pts[outer[k]], m.pts[outer[k2]], yc
        )
        theta_split = (
            _theta_cell_split(hit[3], N_THETA) if hit is not None else None
        )
        add_ogrid_sector_k_baffle(
            b_rout, b_iso, t_rout, t_iso, N_RAD_ISO, nz, k,
            (G_RAD_ISO, 1.0, grading_outer[2]),
            theta_split=theta_split,
        )
        n_layers = len(b_outer) - 1
        ns = _split_radial(N_RAD_COLD, n_layers)
        for i, nr in enumerate(ns):
            add_ogrid_sector_k_baffle(
                b_outer[i], b_outer[i + 1],
                t_outer[i], t_outer[i + 1],
                nr, nz, k, grading_outer, theta_split=theta_split,
            )

    # ---- axial slabs ----
    for iz in range(len(zs) - 1):
        z0, z1 = zs[iz], zs[iz + 1]
        kind = layer_kind(z0, z1)
        baff_i = None
        if isinstance(kind, tuple):
            kind, baff_i = kind
        zm = 0.5 * (z0 + z1)
        nz = nz_for(kind, zm)

        in_shell = kind in ("shell_gap", "baffle")
        in_sheet = kind in ("tube_sheet_in", "tube_sheet_out")
        in_header_cav = kind in ("header_cav_in", "header_cav_out")
        in_header_wall = kind in ("header_wall_in", "header_wall_out")
        in_hot_port = kind in ("hot_port_in", "hot_port_out")

        # ================================================================
        # 1–2. Tube fluid cylinders + hollow tube walls (CAD make_tube_*)
        # ================================================================
        for it, (cx, cy) in enumerate(centers):
            b_core, t_core = tube_core[iz][it], tube_core[iz + 1][it]
            b_rin, t_rin = tube_rin[iz][it], tube_rin[iz + 1][it]
            b_rout, t_rout = tube_rout[iz][it], tube_rout[iz + 1][it]
            b_rhole, t_rhole = tube_rhole[iz][it], tube_rhole[iz + 1][it]
            b_iso, t_iso = tube_iso[iz][it], tube_iso[iz + 1][it]
            b_cell, t_cell = tube_cell[iz][it], tube_cell[iz + 1][it]
            b_hot = [b_core] + [tube_hot_b[iz][ib][it] for ib in range(len(BLEND_TS))] + [b_rin]
            t_hot = [t_core] + [tube_hot_b[iz + 1][ib][it] for ib in range(len(BLEND_TS))] + [t_rin]
            # Cold: OD (or hole) → iso is added separately; outer morph is iso→cell
            cold0_b = b_rout if HOLE_COLLAPSED else b_rhole
            cold0_t = t_rout if HOLE_COLLAPSED else t_rhole
            b_cold_outer = (
                [b_iso]
                + [tube_cold_b[iz][ib][it] for ib in range(len(BLEND_TS))]
                + [b_cell]
            )
            t_cold_outer = (
                [t_iso]
                + [tube_cold_b[iz + 1][ib][it] for ib in range(len(BLEND_TS))]
                + [t_cell]
            )

            g_hot = (G_RAD_HOT, 1.0, 1.0)
            g_cold = (G_RAD_COLD, 1.0, 1.0)
            g_z = (1.0, 1.0, G_Z_BAFFLE) if kind == "baffle" else (1.0, 1.0, 1.0)

            if in_shell or in_sheet:
                # TubeFluid cylinder with smooth square→circle morph
                add_prism("hot_fluid", b_core, t_core, (N_CORE, N_CORE), nz, None, g_z)
                add_morph_chain(
                    "hot_fluid", b_hot, t_hot, N_RAD_HOT, nz, None,
                    (g_hot[0], g_hot[1], g_z[2]),
                )
                # Tube wall hollow cylinder
                add_ogrid_sector(
                    "solid", b_rin, b_rout, t_rin, t_rout, N_WALL, nz, None, g_z
                )

                if in_sheet:
                    if not HOLE_COLLAPSED:
                        add_ogrid_sector(
                            "solid", b_rout, b_rhole, t_rout, t_rhole,
                            N_RAD_HOLE, nz, None, g_cold,
                        )
                    add_cold_ogrid(
                        "solid", cold0_b, cold0_t, b_iso, t_iso,
                        b_cold_outer, t_cold_outer, nz, None,
                        (g_cold[0], g_cold[1], g_z[2]),
                    )
                elif kind == "baffle":
                    if not HOLE_COLLAPSED:
                        add_ogrid_sector(
                            "cold_fluid", b_rout, b_rhole, t_rout, t_rhole,
                            N_RAD_HOLE, nz, None, g_cold,
                        )
                    for k in range(4):
                        add_cold_ogrid_k_baffle(
                            cold0_b, cold0_t, b_iso, t_iso,
                            b_cold_outer, t_cold_outer, nz, k,
                            (g_cold[0], g_cold[1], g_z[2]),
                        )
                else:
                    if not HOLE_COLLAPSED:
                        add_ogrid_sector(
                            "cold_fluid", b_rout, b_rhole, t_rout, t_rhole,
                            N_RAD_HOLE, nz, None, g_cold,
                        )
                    add_cold_ogrid(
                        "cold_fluid", cold0_b, cold0_t, b_iso, t_iso,
                        b_cold_outer, t_cold_outer, nz, None,
                        (g_cold[0], g_cold[1], g_z[2]),
                    )

            elif in_header_cav or in_header_wall or in_hot_port:
                fbot = None
                ftop = None
                if kind == "hot_port_in":
                    fbot = [(
                        "inlet_hot",
                        (b_core[0], b_core[3], b_core[2], b_core[1]),
                    )]
                if kind == "hot_port_out":
                    ftop = [(
                        "outlet_hot",
                        (t_core[0], t_core[1], t_core[2], t_core[3]),
                    )]
                add_prism(
                    "hot_fluid", b_core, t_core, (N_CORE, N_CORE), nz, fbot or ftop
                )
                fmap = None
                if kind == "hot_port_in":
                    fmap = [("inlet_hot", "bot")]
                elif kind == "hot_port_out":
                    fmap = [("outlet_hot", "top")]
                add_morph_chain(
                    "hot_fluid", b_hot, t_hot, N_RAD_HOT, nz, fmap, g_hot
                )
                add_ogrid_sector(
                    "hot_fluid", b_rin, b_rout, t_rin, t_rout, N_WALL, nz
                )
                if not HOLE_COLLAPSED:
                    add_ogrid_sector(
                        "hot_fluid", b_rout, b_rhole, t_rout, t_rhole,
                        N_RAD_HOLE, nz, None, g_cold,
                    )
                add_cold_ogrid(
                    "hot_fluid", cold0_b, cold0_t, b_iso, t_iso,
                    b_cold_outer, t_cold_outer, nz, None, g_cold,
                )

        # ================================================================
        # Bank → shell / header: CONFORMAL collar.
        # Exterior unit-cell verts lie on the pack circle (or, for corner
        # tubes, the outer diagonal is tube-radial). Every exterior
        # pack→shell sector uses a mid-split collar (same wall layers).
        # ================================================================
        def tube_index(it):
            return divmod(it, N_SIDE)  # row, col

        def to_radius(x, y, radius, z):
            r = math.hypot(x, y)
            if r < 1e-14:
                return m.add(radius, 0.0, z)
            return m.add(radius * x / r, radius * y / r, z)

        if in_shell or in_sheet or in_header_cav or in_header_wall or in_hot_port:
            # Outer radius of the fluid fill collar. Hot ports keep the same
            # header radius so they share faces with the header wall (a
            # separate R_NOZZLE O-grid was disconnected and would invert
            # radially because R_NOZZLE < bank half-width).
            if in_header_cav or in_header_wall or in_hot_port:
                r_out_target = r_header
                zone_fill = "hot_fluid"
            elif in_sheet:
                r_out_target = R_SHELL
                zone_fill = "solid"
            else:
                r_out_target = R_SHELL
                zone_fill = None  # decided per sector (baffle window)

            for it, (cx, cy) in enumerate(centers):
                row, col = tube_index(it)
                diag_k = _outer_diagonal_k(row, col)
                # outward sectors: 0=-Y-ish (SW-SE), 1=+X (SE-NE), 2=+Y (NE-NW), 3=-X (NW-SW)
                exterior = []
                if row == 0:
                    exterior.append(0)  # -Y side of pack
                if row == N_SIDE - 1:
                    exterior.append(2)  # +Y
                if col == 0:
                    exterior.append(3)  # -X
                if col == N_SIDE - 1:
                    exterior.append(1)  # +X

                b_cell, t_cell = tube_cell[iz][it], tube_cell[iz + 1][it]
                for k in exterior:
                    k2 = (k + 1) % 4
                    # shell/header points from the same unit-cell corners
                    bx0, by0 = m.pts[b_cell[k]][0], m.pts[b_cell[k]][1]
                    bx1, by1 = m.pts[b_cell[k2]][0], m.pts[b_cell[k2]][1]
                    tx0, ty0 = m.pts[t_cell[k]][0], m.pts[t_cell[k]][1]
                    tx1, ty1 = m.pts[t_cell[k2]][0], m.pts[t_cell[k2]][1]
                    b_a = to_radius(bx0, by0, r_out_target, z0)
                    b_b = to_radius(bx1, by1, r_out_target, z0)
                    t_a = to_radius(tx0, ty0, r_out_target, z1)
                    t_b = to_radius(tx1, ty1, r_out_target, z1)
                    # Keep ray-based shell-ID corners for nozzle wall_fill even
                    # after open_in/out overwrites b_a/b_b with the rect duct.
                    b_a_ray, b_b_ray = b_a, b_b
                    t_a_ray, t_b_ray = t_a, t_b
                    # arcs on outer circle between a-b
                    def _arc_mid(xa, ya, xb, yb, radius, z):
                        a0 = math.atan2(ya, xa)
                        a1 = math.atan2(yb, xb)
                        # shortest directed span
                        da = (a1 - a0 + math.pi) % (2.0 * math.pi) - math.pi
                        am = a0 + 0.5 * da
                        return (radius * math.cos(am), radius * math.sin(am), z)
                    mx, my, mz = _arc_mid(bx0, by0, bx1, by1, r_out_target, z0)
                    m.arc(b_a, b_b, (mx, my, mz))
                    mx, my, mz = _arc_mid(tx0, ty0, tx1, ty1, r_out_target, z1)
                    m.arc(t_a, t_b, (mx, my, mz))

                    zoned = zone_fill
                    if zoned is None:
                        if kind == "baffle":
                            zoned = None  # chord-aware path below
                        else:
                            zoned = "cold_fluid"

                    open_in = (
                        in_shell and kind == "shell_gap"
                        and nozzle_active(zm, "inlet") and k == 2
                        and abs(0.5 * (bx0 + bx1)) < 0.5 * S_L
                    )
                    open_out = (
                        in_shell and kind == "shell_gap"
                        and nozzle_active(zm, "outlet") and k == 0
                        and abs(0.5 * (bx0 + bx1)) < 0.5 * S_L
                    )
                    y_dir_noz = 1.0 if open_in else -1.0

                    def _cyl_noz(x, z, radius, _yd=y_dir_noz):
                        x = max(-(radius - 1e-9), min(radius - 1e-9, x))
                        y = _yd * math.sqrt(max(radius * radius - x * x, 1e-16))
                        return m.add(x, y, z)

                    def _pair_ring(ia, ib, kk):
                        r = [0, 0, 0, 0]
                        r[kk] = ia
                        r[(kk + 1) % 4] = ib
                        return r

                    g_shell = (G_RAD_SHELL, 1.0, G_Z_BAFFLE if kind == "baffle" else 1.0)
                    g_corner = (G_RAD_CORNER, 1.0, G_Z_BAFFLE if kind == "baffle" else 1.0)
                    g_wall = (G_SHELL_WALL, 1.0, G_Z_BAFFLE if kind == "baffle" else 1.0)

                    # Two-layer collar on every exterior pack→shell sector
                    # (same N_RAD_CORNER + N_RAD_SHELL as corner tubes).
                    # Skip nozzles / baffle chord path.
                    split_collar = (
                        not open_in
                        and not open_out
                        and kind != "baffle"
                        and r_out_target == R_SHELL
                    )

                    # Match unit-cell outer face: (N_THETA, nz) on the tube side.
                    # Nozzle: keep FULL pack-edge verts (shared with bank) so the
                    # collar stays face-connected; collapse outer face to the
                    # W×L rectangle at R_SHELL. Wall + stub stay constant W×L.
                    # Rect end-points follow pack-edge X order to avoid inverted hexes.
                    if open_in or open_out:
                        half_x = HALF_X_NOZZLE
                        if m.pts[b_cell[k]][0] <= m.pts[b_cell[k2]][0]:
                            xa, xb = -half_x, +half_x
                        else:
                            xa, xb = +half_x, -half_x
                        b_a = _cyl_noz(xa, z0, r_out_target)
                        b_b = _cyl_noz(xb, z0, r_out_target)
                        t_a = _cyl_noz(xa, z1, r_out_target)
                        t_b = _cyl_noz(xb, z1, r_out_target)
                        v = (
                            b_cell[k], b_a, b_b, b_cell[k2],
                            t_cell[k], t_a, t_b, t_cell[k2],
                        )
                        m.block(zoned, v, (N_RAD_SHELL, N_THETA, nz), None, g_shell)
                    elif kind == "baffle" and zoned is None:
                        # Straight chord cut: split straddling pack→shell sectors.
                        # θ cell counts from pack edge so they match tube cold morph.
                        hit_pack = edge_intersect_chord_y(
                            m.pts[b_cell[k]], m.pts[b_cell[k2]],
                            baffle_chord_y(baff_i),
                        )
                        th_split = (
                            _theta_cell_split(hit_pack[3], N_THETA)
                            if hit_pack is not None else None
                        )
                        add_ogrid_sector_k_baffle(
                            _pair_ring(b_cell[k], b_cell[k2], k),
                            _pair_ring(b_a, b_b, k),
                            _pair_ring(t_cell[k], t_cell[k2], k),
                            _pair_ring(t_a, t_b, k),
                            N_RAD_SHELL, nz, k, g_shell,
                            outer_radius=r_out_target,
                            theta_split=th_split,
                        )
                    elif split_collar:
                        # Two radial layers on each corner ray so neither edge
                        # collapses (shared pack verts stay on r_pack; inset
                        # diagonal has a longer cell→shell span).
                        def _mid_on_ray(pid, z):
                            x, y = m.pts[pid][0], m.pts[pid][1]
                            r0 = math.hypot(x, y)
                            r_m = 0.5 * (r0 + r_out_target)
                            return to_radius(x, y, r_m, z)

                        b_ma = _mid_on_ray(b_cell[k], z0)
                        b_mb = _mid_on_ray(b_cell[k2], z0)
                        t_ma = _mid_on_ray(t_cell[k], z1)
                        t_mb = _mid_on_ray(t_cell[k2], z1)
                        # Outer shell verts stay on the same rays as cell corners
                        b_a = to_radius(bx0, by0, r_out_target, z0)
                        b_b = to_radius(bx1, by1, r_out_target, z0)
                        t_a = to_radius(tx0, ty0, r_out_target, z1)
                        t_b = to_radius(tx1, ty1, r_out_target, z1)
                        mx, my, mz = _arc_mid(bx0, by0, bx1, by1, r_out_target, z0)
                        m.arc(b_a, b_b, (mx, my, mz))
                        mx, my, mz = _arc_mid(tx0, ty0, tx1, ty1, r_out_target, z1)
                        m.arc(t_a, t_b, (mx, my, mz))
                        # Mild arc on the mid ring (lerp of endpoint radii)
                        def _mid_ring_arc(ia, ib, z):
                            xa, ya = m.pts[ia][0], m.pts[ia][1]
                            xb, yb = m.pts[ib][0], m.pts[ib][1]
                            ra, rb = math.hypot(xa, ya), math.hypot(xb, yb)
                            a0 = math.atan2(ya, xa)
                            a1 = math.atan2(yb, xb)
                            da = (a1 - a0 + math.pi) % (2.0 * math.pi) - math.pi
                            am = a0 + 0.5 * da
                            rm = 0.5 * (ra + rb)
                            return (rm * math.cos(am), rm * math.sin(am), z)
                        mx, my, mz = _mid_ring_arc(b_ma, b_mb, z0)
                        m.arc(b_ma, b_mb, (mx, my, mz))
                        mx, my, mz = _mid_ring_arc(t_ma, t_mb, z1)
                        m.arc(t_ma, t_mb, (mx, my, mz))
                        # Layer 1: cell → mid
                        v1 = (
                            b_cell[k], b_ma, b_mb, b_cell[k2],
                            t_cell[k], t_ma, t_mb, t_cell[k2],
                        )
                        m.block(
                            zoned, v1, (N_RAD_CORNER, N_THETA, nz), None, g_corner
                        )
                        # Layer 2: mid → shell
                        v2 = (
                            b_ma, b_a, b_b, b_mb,
                            t_ma, t_a, t_b, t_mb,
                        )
                        m.block(
                            zoned, v2, (N_RAD_SHELL, N_THETA, nz), None, g_shell
                        )
                    else:
                        v = (
                            b_cell[k], b_a, b_b, b_cell[k2],
                            t_cell[k], t_a, t_b, t_cell[k2],
                        )
                        m.block(zoned, v, (N_RAD_SHELL, N_THETA, nz), None, g_shell)

                    # Shell OD wall from the SAME outer verts (conformal)
                    if in_shell or in_sheet:
                        b_oa = to_radius(bx0, by0, R_SHELL_OUT, z0)
                        b_ob = to_radius(bx1, by1, R_SHELL_OUT, z0)
                        t_oa = to_radius(tx0, ty0, R_SHELL_OUT, z1)
                        t_ob = to_radius(tx1, ty1, R_SHELL_OUT, z1)
                        mx, my, mz = _arc_mid(
                            bx0, by0, bx1, by1, R_SHELL_OUT, z0
                        )
                        m.arc(b_oa, b_ob, (mx, my, mz))
                        mx, my, mz = _arc_mid(
                            tx0, ty0, tx1, ty1, R_SHELL_OUT, z1
                        )
                        m.arc(t_oa, t_ob, (mx, my, mz))

                        vw = (
                            b_a, b_oa, b_ob, b_b,
                            t_a, t_oa, t_ob, t_b,
                        )
                        if open_in or open_out:
                            # Constant rectangular duct through wall + stub.
                            # Cross-section W_x × L_z = A_circle along all of Y
                            # until shell entry; N_NOZZLE_Y cells in Y.
                            half_x = HALF_X_NOZZLE
                            y_dir = y_dir_noz
                            # Same X order as collar outer face (xa → xb)
                            bo0 = _cyl_noz(xa, z0, R_SHELL_OUT)
                            bo1 = _cyl_noz(xb, z0, R_SHELL_OUT)
                            to0 = _cyl_noz(xa, z1, R_SHELL_OUT)
                            to1 = _cyl_noz(xb, z1, R_SHELL_OUT)
                            # b_a.. already rect at R_SHELL from collar
                            vw_r = (
                                b_a, bo0, bo1, b_b,
                                t_a, to0, to1, t_b,
                            )
                            m.block(
                                "cold_fluid",
                                vw_r,
                                (N_SHELL_WALL, N_THETA, nz),
                            )

                            # Solid fillers beside nozzle (shell metal + OD).
                            # Outer ends MUST reuse the ray-based shell ID/OD
                            # corners so they share faces with neighbouring
                            # tube shell-wall blocks. (Using _cyl_noz(x_od) on
                            # the ID ring made different points and left 4
                            # disconnected solid islands.)
                            if m.pts[b_oa][0] <= m.pts[b_ob][0]:
                                L_id_b, L_od_b = b_a_ray, b_oa
                                L_id_t, L_od_t = t_a_ray, t_oa
                                R_id_b, R_od_b = b_b_ray, b_ob
                                R_id_t, R_od_t = t_b_ray, t_ob
                            else:
                                L_id_b, L_od_b = b_b_ray, b_ob
                                L_id_t, L_od_t = t_b_ray, t_ob
                                R_id_b, R_od_b = b_a_ray, b_oa
                                R_id_t, R_od_t = t_a_ray, t_oa

                            xL = m.pts[L_od_b][0]
                            xR = m.pts[R_od_b][0]
                            n_span = max(xR - xL, 1e-9)
                            n_left = max(
                                1,
                                int(round(N_THETA * max(-half_x - xL, 0.0) / n_span)),
                            )
                            n_right = max(
                                1,
                                int(round(N_THETA * max(xR - half_x, 0.0) / n_span)),
                            )

                            def _wall_fill_left():
                                if -half_x <= xL + 1e-12 or n_left < 1:
                                    return
                                qi1 = _cyl_noz(-half_x, z0, R_SHELL)
                                qo1 = _cyl_noz(-half_x, z0, R_SHELL_OUT)
                                si1 = _cyl_noz(-half_x, z1, R_SHELL)
                                so1 = _cyl_noz(-half_x, z1, R_SHELL_OUT)
                                m.block(
                                    "solid",
                                    (
                                        L_id_b, L_od_b, qo1, qi1,
                                        L_id_t, L_od_t, so1, si1,
                                    ),
                                    (N_SHELL_WALL, n_left, nz),
                                    [(
                                        "shellOuter",
                                        (L_od_b, qo1, so1, L_od_t),
                                    )],
                                )

                            def _wall_fill_right():
                                if xR <= half_x + 1e-12 or n_right < 1:
                                    return
                                qi0 = _cyl_noz(+half_x, z0, R_SHELL)
                                qo0 = _cyl_noz(+half_x, z0, R_SHELL_OUT)
                                si0 = _cyl_noz(+half_x, z1, R_SHELL)
                                so0 = _cyl_noz(+half_x, z1, R_SHELL_OUT)
                                m.block(
                                    "solid",
                                    (
                                        qi0, qo0, R_od_b, R_id_b,
                                        si0, so0, R_od_t, R_id_t,
                                    ),
                                    (N_SHELL_WALL, n_right, nz),
                                    [(
                                        "shellOuter",
                                        (qo0, R_od_b, R_od_t, so0),
                                    )],
                                )

                            _wall_fill_left()
                            _wall_fill_right()

                            patch = (
                                "inlet_cold" if open_in else "outlet_cold"
                            )
                            face0 = [bo0, bo1, to1, to0]

                            def _shift_y(pid, frac, _yd=y_dir):
                                x, y, z = m.pts[pid]
                                return m.add(
                                    x, y + _yd * NOZZLE_LEN * frac, z
                                )

                            for iy in range(N_NOZZLE_Y):
                                f0 = float(iy) / float(N_NOZZLE_Y)
                                f1 = float(iy + 1) / float(N_NOZZLE_Y)
                                bot = (
                                    face0 if iy == 0
                                    else [_shift_y(pid, f0) for pid in face0]
                                )
                                top = [_shift_y(pid, f1) for pid in face0]
                                if y_dir > 0:
                                    vd = (
                                        bot[0], bot[1], bot[2], bot[3],
                                        top[0], top[1], top[2], top[3],
                                    )
                                else:
                                    vd = (
                                        top[0], top[1], top[2], top[3],
                                        bot[0], bot[1], bot[2], bot[3],
                                    )
                                fmap = None
                                if iy == N_NOZZLE_Y - 1:
                                    fmap = [(
                                        patch,
                                        (top[0], top[1], top[2], top[3]),
                                    )]
                                m.block(
                                    "cold_fluid",
                                    vd,
                                    (N_THETA, nz, 1),
                                    fmap,
                                )

                        elif kind == "baffle":
                            # Split OD with collar for conformal chord faces.
                            # Same θ split as pack edge / collar.
                            hit_pack = edge_intersect_chord_y(
                                m.pts[b_cell[k]], m.pts[b_cell[k2]],
                                baffle_chord_y(baff_i),
                            )
                            th_split = (
                                _theta_cell_split(hit_pack[3], N_THETA)
                                if hit_pack is not None else None
                            )
                            add_ogrid_sector_k_baffle(
                                _pair_ring(b_a, b_b, k),
                                _pair_ring(b_oa, b_ob, k),
                                _pair_ring(t_a, t_b, k),
                                _pair_ring(t_oa, t_ob, k),
                                N_SHELL_WALL, nz, k, g_wall,
                                outer_radius=R_SHELL_OUT,
                                force_zone="solid",
                                face_map_outer="shellOuter",
                                theta_split=th_split,
                            )

                        else:
                            fmap = None
                            if in_shell:
                                fmap = [(
                                    "shellOuter",
                                    (b_oa, b_ob, t_ob, t_oa),
                                )]
                            m.block(
                                "solid",
                                vw,
                                (N_SHELL_WALL, N_THETA, nz),
                                fmap,
                                g_wall,
                            )

                    elif in_header_wall or in_header_cav:
                        # Header cavity + end-wall OD: r_header → shell OD
                        # (cavity must carry the solid shell so end-walls
                        # stay face-connected to the main shell wall)
                        b_oa = to_radius(bx0, by0, R_SHELL_OUT, z0)
                        b_ob = to_radius(bx1, by1, R_SHELL_OUT, z0)
                        t_oa = to_radius(tx0, ty0, R_SHELL_OUT, z1)
                        t_ob = to_radius(tx1, ty1, R_SHELL_OUT, z1)
                        mx, my, mz = _arc_mid(
                            bx0, by0, bx1, by1, R_SHELL_OUT, z0
                        )
                        m.arc(b_oa, b_ob, (mx, my, mz))
                        mx, my, mz = _arc_mid(
                            tx0, ty0, tx1, ty1, R_SHELL_OUT, z1
                        )
                        m.arc(t_oa, t_ob, (mx, my, mz))
                        m.block(
                            "solid",
                            (b_a, b_oa, b_ob, b_b, t_a, t_oa, t_ob, t_b),
                            (N_SHELL_WALL, N_THETA, nz),
                            None,
                            g_wall,
                        )

                    elif in_hot_port:
                        # Port pipe wall around the header-radius bore
                        r_pw = r_header + NOZZLE_WALL
                        b_oa = to_radius(bx0, by0, r_pw, z0)
                        b_ob = to_radius(bx1, by1, r_pw, z0)
                        t_oa = to_radius(tx0, ty0, r_pw, z1)
                        t_ob = to_radius(tx1, ty1, r_pw, z1)
                        mx, my, mz = _arc_mid(bx0, by0, bx1, by1, r_pw, z0)
                        m.arc(b_oa, b_ob, (mx, my, mz))
                        mx, my, mz = _arc_mid(tx0, ty0, tx1, ty1, r_pw, z1)
                        m.arc(t_oa, t_ob, (mx, my, mz))
                        m.block(
                            "solid",
                            (b_a, b_oa, b_ob, b_b, t_a, t_oa, t_ob, t_b),
                            (N_SHELL_WALL, N_THETA, nz),
                            None,
                            g_wall,
                        )

    assert_conformal(m)
    return m, zs, centers, r_header


def assert_conformal(mesh):
    """Fail fast if two blocks share an edge but disagree on its subdivision."""
    # Hex edges → division count (nx along 0-1/3-2, ny along 0-3/1-2, nz along 0-4..)
    edge_pairs_n = (
        ((0, 1), 0), ((3, 2), 0),
        ((0, 3), 1), ((1, 2), 1),
        ((0, 4), 2), ((1, 5), 2), ((2, 6), 2), ((3, 7), 2),
        ((4, 5), 0), ((7, 6), 0),
        ((4, 7), 1), ((5, 6), 1),
    )
    owner = {}
    conflicts = []
    for bi, (_zone, v, div, _gr) in enumerate(mesh.blocks):
        for (a, b), axis in edge_pairs_n:
            i, j = v[a], v[b]
            key = (i, j) if i < j else (j, i)
            n = div[axis]
            if key in owner:
                bj, nj = owner[key]
                if n != nj:
                    conflicts.append((bi, bj, key, n, nj))
            else:
                owner[key] = (bi, n)
    if conflicts:
        msg = ["Non-conformal block interfaces (edge subdivision mismatch):"]
        for bi, bj, key, a, b in conflicts[:20]:
            msg.append(
                "  blocks %d vs %d edge %s: %d vs %d cells" % (bi, bj, key, a, b)
            )
        if len(conflicts) > 20:
            msg.append("  ... and %d more" % (len(conflicts) - 20))
        raise RuntimeError("\n".join(msg))


def write_dict(path, mesh, zs):
    block_edges = set()
    edge_pairs = (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    )
    for _, v, _, _gr in mesh.blocks:
        for a, b in edge_pairs:
            i, j = v[a], v[b]
            block_edges.add((i, j) if i < j else (j, i))
    arcs = {e: mid for e, mid in mesh.arcs.items() if e in block_edges}

    lines = []
    a = lines.append
    a("/*--------------------------------*- C++ -*----------------------------------*\\")
    a("| blockMeshDict — Zahid lab SG-STHX (hot_fluid / cold_fluid / solid)     |")
    a("| Generated by scripts/makeStructuredBlockMesh.py — do not edit by hand.     |")
    a("\\*---------------------------------------------------------------------------*/")
    a("FoamFile")
    a("{")
    a("    version     2.0;")
    a("    format      ascii;")
    a("    class       dictionary;")
    a("    object      blockMeshDict;")
    a("}")
    a("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")
    a("")
    a("scale   1;")
    a("")
    a("vertices")
    a("(")
    for x, y, z in mesh.pts:
        a("    (%.10g %.10g %.10g)" % (x, y, z))
    a(");")
    a("")
    a("blocks")
    a("(")
    for zone, v, div, grading in mesh.blocks:
        a(
            "    hex (%d %d %d %d %d %d %d %d) %s (%d %d %d) "
            "simpleGrading (%.6g %.6g %.6g)"
            % (
                v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], zone,
                div[0], div[1], div[2],
                grading[0], grading[1], grading[2],
            )
        )
    a(");")
    a("")
    a("edges")
    a("(")
    for (i, j), mid in sorted(arcs.items()):
        a("    arc %d %d (%.10g %.10g %.10g)" % (i, j, mid[0], mid[1], mid[2]))
    a(");")
    a("")
    a("boundary")
    a("(")
    for name in ("inlet_hot", "outlet_hot", "inlet_cold", "outlet_cold", "shellOuter"):
        faces = mesh.faces.get(name, [])
        a("    %s" % name)
        a("    {")
        a("        type patch;")
        a("        faces")
        a("        (")
        for f in faces:
            a("            (%d %d %d %d)" % (f[0], f[1], f[2], f[3]))
        a("        );")
        a("    }")
    a(");")
    a("")
    a("defaultPatch")
    a("{")
    a("    type wall;")
    a("    name defaultFaces;")
    a("}")
    a("")
    a("mergePatchPairs")
    a("(")
    a(");")
    a("")
    a("// ************************************************************************* //")
    a("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    case = os.path.dirname(here)
    out = os.path.join(case, "system", "blockMeshDict")

    mesh, zs, centers, r_header = build()
    write_dict(out, mesh, zs)

    counts = {}
    cells = {}
    for zone, _v, div, _gr in mesh.blocks:
        counts[zone] = counts.get(zone, 0) + 1
        cells[zone] = cells.get(zone, 0) + div[0] * div[1] * div[2]

    print("wrote", out)
    print("stations=%d verts=%d blocks=%d headerR=%.4f"
          % (len(zs), len(mesh.pts), len(mesh.blocks), r_header))
    print("  R_HOLE=%.6g R_OUT=%.6g HOLE_COLLAPSED=%s clearance=%.6g"
          % (R_HOLE, R_OUT, HOLE_COLLAPSED, MESH_HOLE_CLEARANCE))
    print("  tube centers (first 4):", centers[:4])
    for z in ("hot_fluid", "cold_fluid", "solid"):
        print("  %s: blocks=%d cells≈%d" % (z, counts.get(z, 0), cells.get(z, 0)))
    print("  total cells≈%d" % sum(cells.values()))
    for name in ("inlet_hot", "outlet_hot", "inlet_cold", "outlet_cold", "shellOuter"):
        print("  patch %s faces=%d" % (name, len(mesh.faces.get(name, []))))
    # bbox check for nozzles
    ys = [p[1] for p in mesh.pts]
    print("  y-range [%.4f, %.4f] (nozzles beyond ±%.3f)"
          % (min(ys), max(ys), R_SHELL))


if __name__ == "__main__":
    main()
