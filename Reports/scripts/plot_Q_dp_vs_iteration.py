#!/usr/bin/env python3
"""Q and shell Δp vs iteration → CSV + SVG under case Reports/.

Q from FO Tout (prefer controlDict surfaceFieldValue.dat):
  Q_hot  = m_dot * Cp * (Tin_h - Tout_h)
  Q_cold = m_dot * Cp * (Tout_c - Tin_c)
  Tin_h=371.47 K, Tin_c=298 K, Cp=4182 (Fluent water)

Δp = faceAvg(p_rgh, inlet_cold) − faceAvg(p_rgh, outlet_cold)
  from reconstructed time dirs and/or processor*/<time>/cold_fluid/p_rgh.
  Optional docker postProcess fallback if field parse fails.

Usage:
  # from case (m027):
  python3 Reports/scripts/plot_Q_dp_vs_iteration.py

  # from CFD-mesh01 (all m0* with FO, or named cases):
  python3 Reports/scripts/plot_Q_dp_vs_iteration.py
  python3 Reports/scripts/plot_Q_dp_vs_iteration.py m027
"""
from __future__ import annotations

import math
import re
import subprocess
import sys
from pathlib import Path

TIN_H = 371.47
TIN_C = 298.0
CP = 4182.0  # Fluent water-liquid

COLOR_Q_HOT = "#c45c26"
COLOR_Q_COLD = "#2f6fed"
COLOR_DP = "#2a7a4b"

# Digitised Fig. 5 experiment (Zahid) — keyed by m_dot
EXP_Q = {
    0.027: 2230.0,
    0.032: 2300.0,
    0.037: 2400.0,
    0.042: 2520.0,
    0.047: 2650.0,
}
EXP_DP = {
    0.027: 2000.0,
    0.032: 2900.0,
    0.037: 3750.0,
    0.042: 4750.0,
    0.047: 5800.0,
}


def _looks_like_case(p: Path) -> bool:
    return (p / "system").is_dir() and (p / "constant").is_dir()


def _looks_like_parent(p: Path) -> bool:
    return any(p.glob("m0*")) and (p / "Reports").is_dir()


def _has_fo(case: Path) -> bool:
    return bool(
        list(case.glob("postProcessing/**/Tout_outlet_hot/**/surfaceFieldValue.dat"))
        or list(case.glob("postProcessing/**/Tout_outlet_hot/**/*.dat"))
    )


def resolve_cases(argv: list[str] | None = None) -> list[Path]:
    """CASE detection: Path(__file__).parents[2], else cwd; parent → m0* with FO."""
    argv = list(argv if argv is not None else sys.argv[1:])
    root = Path(__file__).resolve().parents[2]
    cwd = Path.cwd()

    named = [a for a in argv if not a.startswith("-")]

    def pick_from_parent(parent: Path) -> list[Path]:
        if named:
            out = []
            for n in named:
                c = parent / n if not Path(n).is_absolute() else Path(n)
                if not c.is_dir():
                    c = parent / Path(n).name
                if c.is_dir():
                    out.append(c.resolve())
                else:
                    print(f"WARNING: case not found: {n}", file=sys.stderr)
            return out
        cases = sorted(
            c for c in parent.glob("m0*") if c.is_dir() and _has_fo(c)
        )
        return [c.resolve() for c in cases]

    if _looks_like_case(root):
        return [root]
    if _looks_like_parent(root):
        return pick_from_parent(root)
    if _looks_like_case(cwd):
        return [cwd.resolve()]
    if _looks_like_parent(cwd):
        return pick_from_parent(cwd)
    # Last resort: treat root as case anyway
    return [root]


def m_dot_from_case(case: Path) -> float:
    """m027 → 0.027; else changeDictionaryDict massFlowRate."""
    m = re.fullmatch(r"m(\d+)", case.name)
    if m:
        # keep leading zeros after 'm' so m027 → 0.027 (not 0.27)
        return float(f"0.{m.group(1)}")

    for cand in (
        case / "system" / "cold_fluid" / "changeDictionaryDict",
        case / "system" / "hot_fluid" / "changeDictionaryDict",
        case / "system" / "changeDictionaryDict",
    ):
        if not cand.is_file():
            continue
        text = cand.read_text(errors="ignore")
        mm = re.search(r"massFlowRate\s+([-+eE0-9.]+)\s*;", text)
        if mm:
            return float(mm.group(1))
    raise ValueError(f"cannot determine m_dot for {case}")


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


def extract_from_function_object(case: Path, fo_name: str) -> list[tuple[float, float]]:
    """Parse FO logs → (iter, T_K)."""
    pp = case / "postProcessing"
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


def patch_face_avg_from_field(case: Path, time: float, patch: str) -> float | None:
    """Face-average p_rgh on patch from reconstructed or processor fields."""
    tname = str(int(time)) if float(time).is_integer() else str(time)
    recon = case / tname / "cold_fluid" / "p_rgh"
    if recon.is_file():
        vals = _parse_patch_values(recon.read_text(errors="ignore"), patch)
        if vals:
            return sum(vals) / len(vals)

    vals_all: list[float] = []
    for pr in sorted(case.glob("processor[0-9]*")):
        f = pr / tname / "cold_fluid" / "p_rgh"
        if not f.is_file():
            # try float-named dir (e.g. 1100.0)
            for t, td in _time_dirs(pr):
                if abs(t - time) < 1e-9:
                    f = td / "cold_fluid" / "p_rgh"
                    break
        if not f.is_file():
            continue
        vals = _parse_patch_values(f.read_text(errors="ignore"), patch)
        if vals:
            vals_all.extend(vals)
    if vals_all:
        return sum(vals_all) / len(vals_all)
    return None


def docker_patch_avg_p_rgh(case: Path, time: float) -> float | None:
    """Optional: docker postProcess areaAverage for one time."""
    tname = str(int(time)) if float(time).is_integer() else str(time)
    # parent folder name under Heat-Exchanger (CFD-mesh01)
    mesh_root = case.parent.name
    try:
        out = subprocess.check_output(
            [
                "bash",
                "-lc",
                f"""
                CID=$(docker ps --filter ancestor=opencfd/openfoam-run:latest --format '{{{{.ID}}}}' | head -1)
                [ -n "$CID" ] || exit 2
                docker exec -u "$(id -u):$(id -g)" \\
                  -w /home/openfoam/Validation-study/OpenFOAM/Heat-Exchanger/{mesh_root}/{case.name} \\
                  "$CID" bash -lc '
                  source /usr/lib/openfoam/openfoam2512/etc/bashrc
                  if [ ! -f {tname}/cold_fluid/p_rgh ]; then
                    reconstructPar -region cold_fluid -time {tname} > /dev/null 2>&1 || true
                  fi
                  pin=$(postProcess -func "patchAverage(name=inlet_cold,p_rgh)" -region cold_fluid -time {tname} 2>&1 | awk "/areaAverage\\(inlet_cold\\) of p_rgh/{{print \\$NF}}")
                  pout=$(postProcess -func "patchAverage(name=outlet_cold,p_rgh)" -region cold_fluid -time {tname} 2>&1 | awk "/areaAverage\\(outlet_cold\\) of p_rgh/{{print \\$NF}}")
                  awk -v a="$pin" -v b="$pout" "BEGIN{{printf \\"%.6f\\\\n\\", a-b}}"
                '
                """,
            ],
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=120,
        )
        return float(out.strip().splitlines()[-1])
    except Exception:
        return None


def extract_dp_series(case: Path, *, try_docker: bool = True) -> list[tuple[float, float]]:
    """Δp vs iteration from all available write times (≥1100 preferred, any >0 OK)."""
    times: set[float] = set()
    for t, td in _time_dirs(case):
        if (td / "cold_fluid" / "p_rgh").is_file():
            times.add(t)
    for pr in case.glob("processor[0-9]*"):
        for t, td in _time_dirs(pr):
            if (td / "cold_fluid" / "p_rgh").is_file():
                times.add(t)

    series: list[tuple[float, float]] = []
    for t in sorted(times):
        if t <= 0:
            continue
        pin = patch_face_avg_from_field(case, t, "inlet_cold")
        pout = patch_face_avg_from_field(case, t, "outlet_cold")
        if pin is not None and pout is not None:
            series.append((t, pin - pout))
            continue
        if try_docker:
            dp = docker_patch_avg_p_rgh(case, t)
            if dp is not None and dp == dp:
                series.append((t, dp))
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


def _downsample(
    pts: list[tuple[float, float]], max_pts: int = 1200
) -> list[tuple[float, float]]:
    if len(pts) <= max_pts:
        return pts
    step = max(1, len(pts) // max_pts)
    out = pts[::step]
    if out[-1] != pts[-1]:
        out = out + [pts[-1]]
    return out


def svg_q_plot(
    q_hot: list[tuple[float, float]],
    q_cold: list[tuple[float, float]],
    path: Path,
    *,
    title: str,
    ref_q: float | None,
    width: int = 720,
    height: int = 420,
) -> None:
    if not q_hot and not q_cold:
        return
    all_t = [t for t, _ in q_hot] + [t for t, _ in q_cold]
    all_v = [v for _, v in q_hot] + [v for _, v in q_cold]
    t0, t1 = min(all_t), max(all_t)
    c1 = max(all_v)
    if ref_q is not None:
        c1 = max(c1, ref_q)
    pad = max(50.0, 0.08 * (c1 or 100.0))
    c0, c1 = 0.0, c1 + pad
    if t1 == t0:
        t1 = t0 + 1e-6

    ml, mr, mt, mb = 72, 24, 40, 52
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
        f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#444">Q [W]</text>',
        f'<text x="{ml+pw/2}" y="{height-10}" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#444">Iteration</text>',
    ]

    for tv in nice_ticks(t0, t1, 6):
        xx = x(tv)
        parts.append(
            f'<line x1="{xx}" y1="{mt}" x2="{xx}" y2="{mt+ph}" stroke="#e8e8e8" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{xx}" y="{mt+ph+16}" text-anchor="middle" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#555">'
            f"{tv:g}</text>"
        )

    for cv in nice_ticks(c0, c1, 5):
        yy = y(cv)
        parts.append(
            f'<line x1="{ml}" y1="{yy}" x2="{ml+pw}" y2="{yy}" stroke="#e8e8e8" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{ml-8}" y="{yy+3}" text-anchor="end" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#555">'
            f"{cv:.0f}</text>"
        )

    if ref_q is not None and c0 <= ref_q <= c1:
        yy = y(ref_q)
        parts.append(
            f'<line x1="{ml}" y1="{yy}" x2="{ml+pw}" y2="{yy}" '
            f'stroke="#888" stroke-width="1" stroke-dasharray="4 3"/>'
        )
        parts.append(
            f'<text x="{ml+pw-4}" y="{yy-4}" text-anchor="end" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#666">'
            f"Q_exp={ref_q:.0f} W</text>"
        )

    ly = mt + 8
    if q_hot:
        pts = _downsample(q_hot)
        d = " ".join(
            f"{'M' if j == 0 else 'L'}{x(t):.2f},{y(v):.2f}"
            for j, (t, v) in enumerate(pts)
        )
        parts.append(
            f'<path d="{d}" fill="none" stroke="{COLOR_Q_HOT}" stroke-width="1.6"/>'
        )
        t_last, v_last = q_hot[-1]
        parts.append(
            f'<line x1="{ml+10}" y1="{ly}" x2="{ml+28}" y2="{ly}" '
            f'stroke="{COLOR_Q_HOT}" stroke-width="2"/>'
        )
        parts.append(
            f'<text x="{ml+34}" y="{ly+3}" font-family="Helvetica,Arial,sans-serif" '
            f'font-size="11" fill="#222">'
            f"Q_hot ({v_last:.1f} W @ iter={t_last:g})</text>"
        )
        ly += 16

    if q_cold:
        pts = _downsample(q_cold)
        d = " ".join(
            f"{'M' if j == 0 else 'L'}{x(t):.2f},{y(v):.2f}"
            for j, (t, v) in enumerate(pts)
        )
        parts.append(
            f'<path d="{d}" fill="none" stroke="{COLOR_Q_COLD}" stroke-width="1.6" '
            f'stroke-dasharray="6 3"/>'
        )
        t_last, v_last = q_cold[-1]
        parts.append(
            f'<line x1="{ml+10}" y1="{ly}" x2="{ml+28}" y2="{ly}" '
            f'stroke="{COLOR_Q_COLD}" stroke-width="2" stroke-dasharray="4 2"/>'
        )
        parts.append(
            f'<text x="{ml+34}" y="{ly+3}" font-family="Helvetica,Arial,sans-serif" '
            f'font-size="11" fill="#222">'
            f"Q_cold ({v_last:.1f} W @ iter={t_last:g})</text>"
        )

    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def svg_dp_plot(
    pts: list[tuple[float, float]],
    path: Path,
    *,
    title: str,
    ref_dp: float | None,
    width: int = 720,
    height: int = 420,
) -> None:
    if not pts:
        return
    t0, t1 = pts[0][0], pts[-1][0]
    vals = [v for _, v in pts]
    c1 = max(vals)
    if ref_dp is not None:
        c1 = max(c1, ref_dp)
    pad = max(50.0, 0.08 * (c1 or 100.0))
    c0, c1 = 0.0, c1 + pad
    if t1 == t0:
        t1 = t0 + 1e-6

    ml, mr, mt, mb = 72, 24, 40, 52
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
        f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#444">Δp [Pa]</text>',
        f'<text x="{ml+pw/2}" y="{height-10}" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#444">Iteration</text>',
    ]

    for tv in nice_ticks(t0, t1, 6):
        xx = x(tv)
        parts.append(
            f'<line x1="{xx}" y1="{mt}" x2="{xx}" y2="{mt+ph}" stroke="#e8e8e8" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{xx}" y="{mt+ph+16}" text-anchor="middle" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#555">'
            f"{tv:g}</text>"
        )

    for cv in nice_ticks(c0, c1, 5):
        yy = y(cv)
        parts.append(
            f'<line x1="{ml}" y1="{yy}" x2="{ml+pw}" y2="{yy}" stroke="#e8e8e8" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{ml-8}" y="{yy+3}" text-anchor="end" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#555">'
            f"{cv:.0f}</text>"
        )

    if ref_dp is not None and c0 <= ref_dp <= c1:
        yy = y(ref_dp)
        parts.append(
            f'<line x1="{ml}" y1="{yy}" x2="{ml+pw}" y2="{yy}" '
            f'stroke="#888" stroke-width="1" stroke-dasharray="4 3"/>'
        )
        parts.append(
            f'<text x="{ml+pw-4}" y="{yy-4}" text-anchor="end" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#666">'
            f"Δp_exp={ref_dp:.0f} Pa</text>"
        )

    pts_d = _downsample(pts)
    d = " ".join(
        f"{'M' if j == 0 else 'L'}{x(t):.2f},{y(v):.2f}"
        for j, (t, v) in enumerate(pts_d)
    )
    parts.append(f'<path d="{d}" fill="none" stroke="{COLOR_DP}" stroke-width="1.6"/>')

    # markers for sparse series
    if len(pts) <= 40:
        for t, v in pts:
            parts.append(
                f'<circle cx="{x(t):.2f}" cy="{y(v):.2f}" r="3" fill="{COLOR_DP}"/>'
            )

    t_last, v_last = pts[-1]
    ly = mt + 8
    parts.append(
        f'<line x1="{ml+10}" y1="{ly}" x2="{ml+28}" y2="{ly}" '
        f'stroke="{COLOR_DP}" stroke-width="2"/>'
    )
    parts.append(
        f'<text x="{ml+34}" y="{ly+3}" font-family="Helvetica,Arial,sans-serif" '
        f'font-size="11" fill="#222">'
        f"Δp_shell ({v_last:.1f} Pa @ iter={t_last:g})</text>"
    )

    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def nearest_exp(table: dict[float, float], m_dot: float) -> float | None:
    if m_dot in table:
        return table[m_dot]
    for k, v in table.items():
        if abs(k - m_dot) < 1e-9:
            return v
    return None


def process_case(case: Path) -> None:
    out_data = case / "Reports" / "data"
    out_fig = case / "Reports" / "figures"
    out_data.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)

    m_dot = m_dot_from_case(case)
    ref_q = nearest_exp(EXP_Q, m_dot)
    ref_dp = nearest_exp(EXP_DP, m_dot)

    tho = extract_from_function_object(case, "Tout_outlet_hot")
    tco = extract_from_function_object(case, "Tout_outlet_cold")
    tho = [(t, T) for t, T in tho if t > 0]
    tco = [(t, T) for t, T in tco if t > 0]

    q_hot = [(t, m_dot * CP * (TIN_H - T)) for t, T in tho]
    q_cold = [(t, m_dot * CP * (T - TIN_C)) for t, T in tco]

    q_csv = out_data / "Q_vs_iteration.csv"
    lines = ["iteration,Q_hot_W,Q_cold_W,Tout_hot_K,Tout_cold_K,m_dot_kg_s"]
    # align by iteration (FO usually same times)
    cold_map = {t: (qc, T) for (t, T), (_, qc) in zip(tco, q_cold)}
    if q_hot:
        for (t, T), (_, qh) in zip(tho, q_hot):
            if t in cold_map:
                qc, Tc = cold_map[t]
                lines.append(
                    f"{t:.6g},{qh:.6f},{qc:.6f},{T:.6f},{Tc:.6f},{m_dot:.6g}"
                )
            else:
                lines.append(
                    f"{t:.6g},{qh:.6f},,{T:.6f},,{m_dot:.6g}"
                )
        # any cold-only times
        hot_times = {t for t, _ in tho}
        for t, T in tco:
            if t not in hot_times:
                qc = m_dot * CP * (T - TIN_C)
                lines.append(f"{t:.6g},,{qc:.6f},,{T:.6f},{m_dot:.6g}")
    elif q_cold:
        for (t, T), (_, qc) in zip(tco, q_cold):
            lines.append(f"{t:.6g},,{qc:.6f},,{T:.6f},{m_dot:.6g}")
    q_csv.write_text("\n".join(lines) + "\n")

    q_svg = out_fig / "Q_vs_iteration.svg"
    if q_hot or q_cold:
        svg_q_plot(
            q_hot,
            q_cold,
            q_svg,
            title=f"Heat transfer rate vs iteration — {case.name} (m={m_dot:g} kg/s)",
            ref_q=ref_q,
        )

    dp = extract_dp_series(case)
    dp_csv = out_data / "dP_vs_iteration.csv"
    dp_lines = ["iteration,dP_shell_Pa,note"]
    for t, v in dp:
        dp_lines.append(f"{t:.6g},{v:.6f},p_rgh(inlet_cold)-p_rgh(outlet_cold)")
    dp_csv.write_text("\n".join(dp_lines) + "\n")

    dp_svg = out_fig / "dP_vs_iteration.svg"
    if dp:
        svg_dp_plot(
            dp,
            dp_svg,
            title=f"Shell Δp vs iteration — {case.name} (m={m_dot:g} kg/s)",
            ref_dp=ref_dp,
        )

    # summary
    def rel(p: Path) -> str:
        try:
            return str(p.relative_to(case))
        except ValueError:
            return str(p)

    if q_hot:
        print(
            f"{case.name}: wrote {rel(q_svg)}  "
            f"({len(q_hot)} samples, last Q_hot={q_hot[-1][1]:.1f} W @ iter={q_hot[-1][0]:g})"
        )
        if q_cold:
            print(
                f"         last Q_cold={q_cold[-1][1]:.1f} W @ iter={q_cold[-1][0]:g}"
            )
    else:
        print(f"{case.name}: WARNING: no Q data (missing FO Tout)")
    print(f"         wrote {rel(q_csv)}")

    if dp:
        print(
            f"{case.name}: wrote {rel(dp_svg)}  "
            f"({len(dp)} samples, last Δp={dp[-1][1]:.1f} Pa @ iter={dp[-1][0]:g})"
        )
    else:
        print(f"{case.name}: WARNING: no Δp data (no p_rgh write times)")
    print(f"         wrote {rel(dp_csv)}")


def main() -> int:
    cases = resolve_cases()
    if not cases:
        print("No cases to process", file=sys.stderr)
        return 1
    for case in cases:
        process_case(case)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
