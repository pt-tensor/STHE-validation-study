#!/usr/bin/env python3
"""Extract outlet_hot / outlet_cold T vs iteration and write SVGs.

Preference order:
  1. postProcessing/**/Tout_outlet_*/surfaceFieldValue.dat (controlDict FO)
  2. reconstructed time dirs / processor*/<time>/<region>/T

  Reports/data/Tout_hot_vs_iteration.csv
  Reports/data/Tout_cold_vs_iteration.csv
  Reports/figures/Tout_hot_vs_iteration.svg
  Reports/figures/Tout_cold_vs_iteration.svg

Usage (from CFD/ or anywhere):
  python3 Reports/scripts/plot_outlet_temperatures.py
"""
from __future__ import annotations

import math
import re
from pathlib import Path

CASE = Path(__file__).resolve().parents[2]  # .../Heat-Exchanger/CFD
OUT_DATA = CASE / "Reports" / "data"
OUT_FIG = CASE / "Reports" / "figures"

# Zahid lab SG: Tin from Fig. 6 legend max/min as inlets (371.47 / 298 K).
# No closed-form Tout benchmark — omit ref lines (None).
REF_TH_OUT_C = None
REF_TC_OUT_C = None

COLOR_HOT = "#c45c26"
COLOR_COLD = "#2f6fed"

# key → (region, patch, ref_C, color, legend, FO name in controlDict)
PATCH = {
    "hot": (
        "hot_fluid",
        "outlet_hot",
        REF_TH_OUT_C,
        COLOR_HOT,
        "T_out,hot",
        "Tout_outlet_hot",
    ),
    "cold": (
        "cold_fluid",
        "outlet_cold",
        REF_TC_OUT_C,
        COLOR_COLD,
        "T_out,cold",
        "Tout_outlet_cold",
    ),
}


def _parse_patch_values(text: str, patch: str) -> list[float]:
    """Return face values for a named patch (empty if patch has 0 faces)."""
    m = re.search(
        rf"{re.escape(patch)}\s*\{{(.*?)\n\s*\}}",
        text,
        flags=re.S,
    )
    if not m:
        return []
    block = m.group(1)
    mu = re.search(r"value\s+uniform\s+([-+eE0-9.]+)\s*;", block)
    if mu:
        return [float(mu.group(1))]
    mn = re.search(
        r"value\s+nonuniform\s+List<scalar>\s*\n\s*(\d+)\s*\n\s*\((.*?)\)\s*;",
        block,
        flags=re.S,
    )
    if not mn:
        if re.search(r"value\s+nonuniform\s+List<scalar>\s+0\(\)\s*;", block):
            return []
        return []
    n = int(mn.group(1))
    body = mn.group(2)
    vals = [float(x) for x in body.split()]
    if n and len(vals) != n:
        pass
    return vals


def _time_dirs(root: Path) -> list[tuple[float, Path]]:
    out: list[tuple[float, Path]] = []
    if not root.exists():
        return out
    for p in root.iterdir():
        if not p.is_dir():
            continue
        try:
            t = float(p.name)
        except ValueError:
            continue
        out.append((t, p))
    out.sort(key=lambda x: x[0])
    return out


def extract_from_function_object(fo_name: str) -> list[tuple[float, float]]:
    """Parse FO logs → (iter, T_K).

    ESI multi-region writes:
      postProcessing/<region>/<fo_name>/0/surfaceFieldValue.dat
    """
    pp = CASE / "postProcessing"
    if not pp.is_dir():
        return []
    dats = sorted(pp.glob(f"**/{fo_name}/**/surfaceFieldValue.dat"))
    if not dats:
        dats = sorted(pp.glob(f"**/{fo_name}/**/*.dat"))
    series: list[tuple[float, float]] = []
    seen: set[float] = set()
    for dat in dats:
        for line in dat.read_text(errors="ignore").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                t = float(parts[0])
                T = float(parts[-1])
            except ValueError:
                continue
            if t in seen:
                continue
            seen.add(t)
            series.append((t, T))
    series.sort(key=lambda x: x[0])
    return series


def extract_series(
    region: str, patch: str, fo_name: str = ""
) -> list[tuple[float, float]]:
    """Return (iteration, T_K); prefer controlDict FO log, else field dumps."""
    if fo_name:
        fo = extract_from_function_object(fo_name)
        if fo:
            return fo

    recon: dict[float, list[float]] = {}
    for t, td in _time_dirs(CASE):
        f = td / region / "T"
        if not f.exists():
            continue
        vals = _parse_patch_values(f.read_text(errors="ignore"), patch)
        if vals:
            recon[t] = vals

    proc_roots = sorted(CASE.glob("processor[0-9]*"))
    proc_times: dict[float, list[float]] = {}
    for pr in proc_roots:
        for t, td in _time_dirs(pr):
            f = td / region / "T"
            if not f.exists():
                continue
            vals = _parse_patch_values(f.read_text(errors="ignore"), patch)
            if not vals:
                continue
            proc_times.setdefault(t, []).extend(vals)

    times = sorted(set(recon) | set(proc_times))
    series: list[tuple[float, float]] = []
    for t in times:
        vals = recon.get(t) or proc_times.get(t) or []
        if not vals:
            continue
        series.append((t, sum(vals) / len(vals)))
    return series


def nice_ticks(a: float, b: float, n: int = 5) -> list[float]:
    span = b - a
    if span <= 0:
        return [a]
    step = 10 ** math.floor(math.log10(span / max(n, 1)))
    for m in (1, 2, 5, 10):
        if span / (m * step) <= n + 1:
            step = m * step
            break
    start = math.ceil(a / step - 1e-12) * step
    ticks: list[float] = []
    v = start
    while v <= b + 1e-9 * max(1.0, abs(b)):
        ticks.append(v)
        v += step
        if len(ticks) > 20:
            break
    return ticks or [a, b]


def svg_plot(
    pts_c: list[tuple[float, float]],
    path: Path,
    *,
    title: str,
    ylabel: str,
    color: str,
    ref_c: float | None,
    ref_label: str,
    series_name: str,
    width: int = 720,
    height: int = 420,
) -> None:
    if not pts_c:
        return

    t0, t1 = pts_c[0][0], pts_c[-1][0]
    vals = [c for _, c in pts_c]
    c1 = max(vals)
    if ref_c is not None:
        c1 = max(c1, ref_c)
    pad = max(1.0, 0.08 * (c1 or 10.0))
    c0, c1 = 0.0, c1 + pad
    if t1 == t0:
        t1 = t0 + 1e-6

    ml, mr, mt, mb = 64, 24, 40, 52
    pw, ph = width - ml - mr, height - mt - mb

    def x(t: float) -> float:
        return ml + (t - t0) / (t1 - t0) * pw

    def y(c: float) -> float:
        c_clamped = min(max(c, c0), c1)
        return mt + (1 - (c_clamped - c0) / (c1 - c0)) * ph

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="24" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="14" fill="#111">{title}</text>',
        f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#333" stroke-width="1"/>',
        f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#333" stroke-width="1"/>',
        f'<text x="{ml-10}" y="{mt-8}" text-anchor="end" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#444">{ylabel}</text>',
        f'<text x="{ml+pw/2}" y="{height-10}" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#444">Iteration</text>',
    ]

    for tv in nice_ticks(t0, t1, 6):
        xx = x(tv)
        parts.append(
            f'<line x1="{xx}" y1="{mt}" x2="{xx}" y2="{mt+ph}" stroke="#e8e8e8" stroke-width="1"/>'
        )
        lab = f"{tv:g}"
        parts.append(
            f'<text x="{xx}" y="{mt+ph+16}" text-anchor="middle" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#555">'
            f"{lab}</text>"
        )

    for cv in nice_ticks(c0, c1, 5):
        yy = y(cv)
        parts.append(
            f'<line x1="{ml}" y1="{yy}" x2="{ml+pw}" y2="{yy}" stroke="#e8e8e8" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{ml-8}" y="{yy+3}" text-anchor="end" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#555">'
            f"{cv:.1f}</text>"
        )

    if ref_c is not None and c0 <= ref_c <= c1:
        yy = y(ref_c)
        parts.append(
            f'<line x1="{ml}" y1="{yy}" x2="{ml+pw}" y2="{yy}" '
            f'stroke="#888" stroke-width="1" stroke-dasharray="4 3"/>'
        )
        parts.append(
            f'<text x="{ml+pw-4}" y="{yy-4}" text-anchor="end" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#666">'
            f"{ref_label}={ref_c:.1f} °C</text>"
        )

    step = max(1, len(pts_c) // 1200)
    pts_d = pts_c[::step]
    if pts_d[-1] != pts_c[-1]:
        pts_d = pts_d + [pts_c[-1]]
    d = " ".join(
        f"{'M' if j == 0 else 'L'}{x(t):.2f},{y(c):.2f}"
        for j, (t, c) in enumerate(pts_d)
    )
    parts.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="1.6"/>')

    t_last, c_last = pts_c[-1]
    ly = mt + 8
    parts.append(
        f'<line x1="{ml+10}" y1="{ly}" x2="{ml+28}" y2="{ly}" '
        f'stroke="{color}" stroke-width="2"/>'
    )
    parts.append(
        f'<text x="{ml+34}" y="{ly+3}" font-family="Helvetica,Arial,sans-serif" '
        f'font-size="11" fill="#222">'
        f"{series_name} ({c_last:.2f} °C @ iter={t_last:g})</text>"
    )

    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def write_csv(path: Path, pts_c: list[tuple[float, float]]) -> None:
    lines = ["iteration,T_out_C,T_out_K"]
    for t, tc in pts_c:
        lines.append(f"{t:.6g},{tc:.6f},{tc + 273.15:.6f}")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)

    for key, (region, patch, ref_c, color, name, fo_name) in PATCH.items():
        series_k = extract_series(region, patch, fo_name)
        if not series_k:
            print(f"WARNING: no data for {patch} ({region})")
            continue
        # Skip pure IC (iter 0) so the curve tracks the steady run
        series_k = [(it, T) for it, T in series_k if it > 0]
        if not series_k:
            print(f"WARNING: no post-IC data for {patch} ({region})")
            continue
        series_c = [(t, T - 273.15) for t, T in series_k]
        csv_path = OUT_DATA / f"Tout_{key}_vs_iteration.csv"
        svg_path = OUT_FIG / f"Tout_{key}_vs_iteration.svg"
        write_csv(csv_path, series_c)
        svg_plot(
            series_c,
            svg_path,
            title=f"Outlet temperature vs iteration — {patch} (CFD)",
            ylabel="T [°C]",
            color=color,
            ref_c=ref_c,
            ref_label="ref",
            series_name=name,
        )
        print(
            f"wrote {svg_path.relative_to(CASE)}  "
            f"({len(series_c)} samples, "
            f"last {series_c[-1][1]:.2f} °C @ iter={series_c[-1][0]:g})"
        )
        print(f"wrote {csv_path.relative_to(CASE)}")


if __name__ == "__main__":
    main()
