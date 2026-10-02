#!/usr/bin/env python3
"""Build Zahid Fig. 5 validation table + SVGs (experiment + OpenFOAM).

Primary reference: experimental points digitised into
  Reference/(MAIN REFERENCE) ….md  → Fig. 5(a)/(b)

Validation focus: experiment vs OpenFOAM only.

Usage:
  python3 Reports/scripts/build_validation_fig5.py
"""
from __future__ import annotations

import csv
import math
import os
import re
from pathlib import Path

CASE = Path(__file__).resolve().parents[2]
OUT_DATA = CASE / "Reports" / "data"
OUT_FIG = CASE / "Reports" / "figures"
OUT_MD = CASE / "Reports" / "validation_fig5.md"
REF_NAME = (
    "(MAIN REFERENCE) Experimental and CFD Simulation Study of Shell and Tube "
    "Heat Exchangers with Different Baffle Segment Configurations.md"
)


def find_hx_root(start: Path) -> Path:
    """Walk up from CASE until Reference/(MAIN REFERENCE)….md exists."""
    for p in [start, *start.parents]:
        cand = p / "Reference" / REF_NAME
        if cand.is_file():
            return p
    return start.parent


HX = find_hx_root(CASE)
REF_MD = HX / "Reference" / REF_NAME


def ref_md_rel_from_reports() -> str:
    """Markdown-relative link from CASE/Reports/ to MAIN REFERENCE."""
    try:
        rel = Path(os.path.relpath(REF_MD, CASE / "Reports")).as_posix()
    except Exception:
        rel = "../../Reference/" + REF_NAME
    return rel.replace(" ", "%20")

# ---------------------------------------------------------------------------
# Table 1 geometry + Tin from Fig. 6 legend (assumed inlets)
# ---------------------------------------------------------------------------
DO = 0.006
DI = 0.004
DS = 0.044
PT = 0.012
B = 0.0356
L = 0.186
NT = 9
NB = 4
KW = 15.2

TIN_H = 371.47  # Fig. 6 legend max [K]
TIN_C = 298.0   # Fig. 6 legend min [K]
M_DOTS = (0.027, 0.032, 0.037, 0.042, 0.047)

PROP_H = dict(rho=960.0, mu=2.9e-4, k=0.677, Cp=4210.0, Pr=1.8)
PROP_C = dict(rho=997.0, mu=8.9e-4, k=0.607, Cp=4182.0, Pr=6.1)

# Fallback if REF_MD parse fails (same numbers as digitised Fig. 5)
EXP_Q_FALLBACK = {
    0.027: 2230.0,
    0.032: 2300.0,
    0.037: 2400.0,
    0.042: 2520.0,
    0.047: 2650.0,
}
EXP_DP_FALLBACK = {
    0.027: 2000.0,
    0.032: 2900.0,
    0.037: 3750.0,
    0.042: 4750.0,
    0.047: 5800.0,
}

CFD_DEFAULT = [
    (
        0.047,
        1260.4,
        238.8,
        53.587,
        31.413,
        "CFD counterflow t=2000; FO Tout; dP=p_rgh",
    ),
]


def parse_fig5_experiment(md_path: Path) -> tuple[dict[float, float], dict[float, float]]:
    """Parse digitised Fig. 5 tables from MAIN REFERENCE md."""
    if not md_path.is_file():
        return dict(EXP_Q_FALLBACK), dict(EXP_DP_FALLBACK)
    text = md_path.read_text(encoding="utf-8")
    # Find Fig. 5 section
    m = re.search(r"Fig\.\s*5\..*?(?=## |\Z)", text, re.S | re.I)
    block = m.group(0) if m else text
    q: dict[float, float] = {}
    dp: dict[float, float] = {}
    # rows like: | 0.047 | 2650 |
    rows = re.findall(
        r"\|\s*(0\.0[2-4][0-9])\s*\|\s*([0-9]+(?:\.[0-9]+)?)\s*\|",
        block,
    )
    # First 5 → Q, next 5 → dP (as edited in the md)
    vals = [(float(a), float(b)) for a, b in rows]
    if len(vals) >= 10:
        for mdot, v in vals[:5]:
            q[mdot] = v
        for mdot, v in vals[5:10]:
            dp[mdot] = v
    elif len(vals) >= 5:
        # only one table found — assume Q
        for mdot, v in vals[:5]:
            q[mdot] = v
        dp = dict(EXP_DP_FALLBACK)
    else:
        return dict(EXP_Q_FALLBACK), dict(EXP_DP_FALLBACK)
    if not q:
        q = dict(EXP_Q_FALLBACK)
    if not dp:
        dp = dict(EXP_DP_FALLBACK)
    return q, dp


def shell_flow_area() -> float:
    return DS * (PT - DO) * B / PT


def shell_De_square() -> float:
    return 4.0 * (PT * PT - math.pi * DO * DO / 4.0) / (math.pi * DO)


def kern_hs(m_dot: float) -> tuple[float, float, float]:
    As = shell_flow_area()
    De = shell_De_square()
    Gs = m_dot / As
    Re = De * Gs / PROP_C["mu"]
    Nu = 0.36 * (Re ** 0.55) * (PROP_C["Pr"] ** (1.0 / 3.0))
    hs = Nu * PROP_C["k"] / De
    return hs, Re, Gs


def tube_hi(m_dot: float) -> tuple[float, float]:
    At = NT * math.pi * DI * DI / 4.0
    Gt = m_dot / At
    Re = DI * Gt / PROP_H["mu"]
    Nu = 0.023 * (Re ** 0.8) * (PROP_H["Pr"] ** 0.3)
    hi = Nu * PROP_H["k"] / DI
    return hi, Re


def overall_Uo(hi: float, ho: float) -> float:
    wall = DO * math.log(DO / DI) / (2.0 * KW)
    return 1.0 / (DO / (DI * hi) + wall + 1.0 / ho)


def heat_transfer_area() -> float:
    return NT * math.pi * DO * L


def kern_Q(m_dot: float) -> dict:
    hs, Re_s, Gs = kern_hs(m_dot)
    hi, Re_t = tube_hi(m_dot)
    Uo = overall_Uo(hi, hs)
    A = heat_transfer_area()
    Ch = m_dot * PROP_H["Cp"]
    Cc = m_dot * PROP_C["Cp"]
    Cmin = min(Ch, Cc)
    Cmax = max(Ch, Cc)
    Cr = Cmin / Cmax
    NTU = Uo * A / Cmin
    if abs(Cr - 1.0) < 1e-9:
        eps = NTU / (1.0 + NTU)
    else:
        eps = (1.0 - math.exp(-NTU * (1.0 - Cr))) / (
            1.0 - Cr * math.exp(-NTU * (1.0 - Cr))
        )
    Qmax = Cmin * (TIN_H - TIN_C)
    Q = eps * Qmax
    return {
        "Q_W": Q,
        "Uo": Uo,
        "hs": hs,
        "hi": hi,
        "Re_s": Re_s,
        "Re_t": Re_t,
        "Gs": Gs,
        "NTU": NTU,
        "eps": eps,
        "Tho_C": (TIN_H - Q / Ch) - 273.15,
        "Tco_C": (TIN_C + Q / Cc) - 273.15,
    }


def esso_dP(m_dot: float) -> dict:
    As = shell_flow_area()
    De = shell_De_square()
    Gs = m_dot / As
    Re = max(De * Gs / PROP_C["mu"], 1.0)
    j_f = math.exp(0.576 - 0.19 * math.log(Re))
    dP = 8.0 * j_f * (DS / De) * (Gs * Gs / (2.0 * PROP_C["rho"])) * (L / B)
    return {"dP_Pa": dP, "f": j_f, "Re_s": Re, "Gs": Gs}


def load_cfd_points() -> list[dict]:
    path = OUT_DATA / "cfd_validation_points.csv"
    rows: list[dict] = []
    if path.is_file():
        with path.open(newline="") as fh:
            for r in csv.DictReader(fh):
                rows.append(
                    {
                        "m_dot": float(r["m_dot_kg_s"]),
                        "Q_W": float(r["Q_W"]),
                        "dP_shell_Pa": float(r["dP_shell_Pa"]),
                        "Tho_C": float(r["Tho_C"]) if r.get("Tho_C") else None,
                        "Tco_C": float(r["Tco_C"]) if r.get("Tco_C") else None,
                        "note": r.get("note", ""),
                    }
                )
        return rows
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["m_dot_kg_s", "Q_W", "dP_shell_Pa", "Tho_C", "Tco_C", "note"])
        for m, Q, dP, Tho, Tco, note in CFD_DEFAULT:
            w.writerow([m, f"{Q:.1f}", f"{dP:.1f}", f"{Tho:.3f}", f"{Tco:.3f}", note])
            rows.append(
                {
                    "m_dot": m,
                    "Q_W": Q,
                    "dP_shell_Pa": dP,
                    "Tho_C": Tho,
                    "Tco_C": Tco,
                    "note": note,
                }
            )
    return rows


def build_table(
    cfd: list[dict],
    exp_q: dict[float, float],
    exp_dp: dict[float, float],
) -> list[dict]:
    cfd_by_m = {round(r["m_dot"], 5): r for r in cfd}
    rows = []
    for m in M_DOTS:
        k = kern_Q(m)
        e = esso_dP(m)
        c = cfd_by_m.get(round(m, 5))
        Q_exp = exp_q.get(m)
        dP_exp = exp_dp.get(m)
        Q_cfd = c["Q_W"] if c else None
        dP_cfd = c["dP_shell_Pa"] if c else None
        err_Q_exp = (
            100.0 * (Q_cfd - Q_exp) / Q_exp
            if Q_cfd is not None and Q_exp
            else None
        )
        err_dP_exp = (
            100.0 * (dP_cfd - dP_exp) / dP_exp
            if dP_cfd is not None and dP_exp
            else None
        )
        err_Q_kern = (
            100.0 * (Q_cfd - k["Q_W"]) / k["Q_W"]
            if Q_cfd is not None and k["Q_W"]
            else None
        )
        err_dP_esso = (
            100.0 * (dP_cfd - e["dP_Pa"]) / e["dP_Pa"]
            if dP_cfd is not None and e["dP_Pa"]
            else None
        )
        rows.append(
            {
                "m_dot": m,
                "Q_exp_W": Q_exp,
                "Q_Kern_W": k["Q_W"],
                "Q_CFD_W": Q_cfd,
                "err_Q_vs_exp_pct": err_Q_exp,
                "err_Q_vs_Kern_pct": err_Q_kern,
                "dP_exp_Pa": dP_exp,
                "dP_Esso_Pa": e["dP_Pa"],
                "dP_CFD_Pa": dP_cfd,
                "err_dP_vs_exp_pct": err_dP_exp,
                "err_dP_vs_Esso_pct": err_dP_esso,
                "hs_Kern": k["hs"],
                "Re_s": k["Re_s"],
                "Tho_CFD_C": c["Tho_C"] if c else None,
                "Tco_CFD_C": c["Tco_C"] if c else None,
            }
        )
    return rows


def write_paper_csv(exp_q: dict[float, float], exp_dp: dict[float, float]) -> Path:
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    path = OUT_DATA / "paper_fig5_experiment.csv"
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["m_dot_kg_s", "Q_exp_W", "dP_exp_Pa", "source"])
        for m in M_DOTS:
            w.writerow(
                [
                    f"{m:.3f}",
                    f"{exp_q[m]:.0f}",
                    f"{exp_dp[m]:.0f}",
                    "MAIN REFERENCE Fig.5 digitised",
                ]
            )
    return path


def write_csv(rows: list[dict]) -> Path:
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    path = OUT_DATA / "validation_fig5.csv"
    cols = [
        "m_dot_kg_s",
        "Q_exp_W",
        "Q_OpenFOAM_W",
        "err_Q_vs_exp_pct",
        "Q_Kern_rebuild_W",
        "err_Q_vs_Kern_pct",
        "dP_exp_Pa",
        "dP_OpenFOAM_Pa",
        "err_dP_vs_exp_pct",
        "dP_Esso_rebuild_Pa",
        "err_dP_vs_Esso_pct",
        "Tho_OpenFOAM_C",
        "Tco_OpenFOAM_C",
    ]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow(
                [
                    f"{r['m_dot']:.3f}",
                    "" if r["Q_exp_W"] is None else f"{r['Q_exp_W']:.0f}",
                    "" if r["Q_CFD_W"] is None else f"{r['Q_CFD_W']:.1f}",
                    ""
                    if r["err_Q_vs_exp_pct"] is None
                    else f"{r['err_Q_vs_exp_pct']:.1f}",
                    f"{r['Q_Kern_W']:.1f}",
                    ""
                    if r["err_Q_vs_Kern_pct"] is None
                    else f"{r['err_Q_vs_Kern_pct']:.1f}",
                    "" if r["dP_exp_Pa"] is None else f"{r['dP_exp_Pa']:.0f}",
                    "" if r["dP_CFD_Pa"] is None else f"{r['dP_CFD_Pa']:.1f}",
                    ""
                    if r["err_dP_vs_exp_pct"] is None
                    else f"{r['err_dP_vs_exp_pct']:.1f}",
                    f"{r['dP_Esso_Pa']:.1f}",
                    ""
                    if r["err_dP_vs_Esso_pct"] is None
                    else f"{r['err_dP_vs_Esso_pct']:.1f}",
                    "" if r["Tho_CFD_C"] is None else f"{r['Tho_CFD_C']:.2f}",
                    "" if r["Tco_CFD_C"] is None else f"{r['Tco_CFD_C']:.2f}",
                ]
            )
    return path


def write_md(rows: list[dict]) -> Path:
    cfd_ms = [r["m_dot"] for r in rows if r["Q_CFD_W"] is not None]
    if len(cfd_ms) == 1:
        m_desc = f"ṁ = {cfd_ms[0]:.3f} kg/s"
    elif cfd_ms:
        m_desc = f"ṁ = {cfd_ms[0]:.3f}–{cfd_ms[-1]:.3f} kg/s ({len(cfd_ms)} cases)"
    else:
        m_desc = "ṁ = (no OF points yet)"
    case_label = CASE.name if CASE.name.startswith("m0") else "CFD/"
    ref_href = ref_md_rel_from_reports()
    lines = [
        "# Model validation — Zahid Fig. 5 (SG-STHX)",
        "",
        "Primary reference: **experimental** points digitised in "
        f"[MAIN REFERENCE]({ref_href}) Fig. 5.",
        "",
        f"**OpenFOAM:** `{case_label}` counterflow, kOmegaSST, {m_desc}, "
        "Tin = 371.47 / 298 K (Fig. 6 legend max/min as inlets).",
        "",
        "Validation focus: **experiment vs OpenFOAM** only.",
        "",
        "## Validation vs experiment",
        "",
        "| ṁ [kg/s] | Q<sub>exp</sub> [W] | Q<sub>OF</sub> [W] | err Q [%] | "
        "Δp<sub>exp</sub> [Pa] | Δp<sub>OF</sub> [Pa] | err Δp [%] |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        q_of = "—" if r["Q_CFD_W"] is None else f"{r['Q_CFD_W']:.1f}"
        e_q = (
            "—"
            if r["err_Q_vs_exp_pct"] is None
            else f"{r['err_Q_vs_exp_pct']:+.1f}"
        )
        dp_of = "—" if r["dP_CFD_Pa"] is None else f"{r['dP_CFD_Pa']:.1f}"
        e_dp = (
            "—"
            if r["err_dP_vs_exp_pct"] is None
            else f"{r['err_dP_vs_exp_pct']:+.1f}"
        )
        lines.append(
            f"| {r['m_dot']:.3f} | {r['Q_exp_W']:.0f} | {q_of} | {e_q} | "
            f"{r['dP_exp_Pa']:.0f} | {dp_of} | {e_dp} |"
        )
    lines += [
        "",
        "## Figures",
        "",
        "- [`figures/validation_heat_transfer_rate.svg`](figures/validation_heat_transfer_rate.svg)",
        "- [`figures/validation_pressure_drop.svg`](figures/validation_pressure_drop.svg)",
        "",
        "## Notes",
        "",
        "- Errors are **OpenFOAM − experiment** / experiment.",
        "- Tin from Fig. 6 colour-bar extremes: **371.47 K** (hot), **298 K** (cold).",
        "- Shell nozzles: equivalent Ø **5.5 mm** (D_shell/8), uniform "
        "W×L stub into the shell (multi-cell in Y).",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    return OUT_MD


def _nice_ticks(a: float, b: float, n: int = 5) -> list[float]:
    if b <= a:
        return [a]
    span = b - a
    raw = span / max(n - 1, 1)
    exp = math.floor(math.log10(raw)) if raw > 0 else 0
    base = 10**exp
    step = base
    for mult in (1, 2, 2.5, 5, 10):
        if mult * base >= raw * 0.8:
            step = mult * base
            break
    start = math.ceil(a / step) * step
    ticks = []
    x = start
    while x <= b + 0.5 * step:
        if a - 1e-12 <= x <= b + 1e-12:
            ticks.append(x)
        x += step
    return ticks or [a, b]


def svg_xy(
    series: list[tuple[str, list[tuple[float, float]], str, str]],
    path: Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    width: int = 720,
    height: int = 440,
) -> None:
    pts_all = [p for _, pts, _, _ in series for p in pts]
    if not pts_all:
        return
    xs = [p[0] for p in pts_all]
    ys = [p[1] for p in pts_all]
    x0, x1 = min(xs), max(xs)
    y1 = max(ys)
    pad_x = 0.05 * (x1 - x0 or 1.0)
    pad_y = 0.10 * (y1 or 1.0)
    x0, x1 = x0 - pad_x, x1 + pad_x
    y0, y1 = 0.0, y1 + pad_y

    ml, mr, mt, mb = 72, 28, 44, 56
    pw, ph = width - ml - mr, height - mt - mb

    def X(x: float) -> float:
        return ml + (x - x0) / (x1 - x0) * pw

    def Y(y: float) -> float:
        return mt + (1.0 - (y - y0) / (y1 - y0)) * ph

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-label="{title}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="26" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="14" fill="#111">{title}</text>',
        f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#333" stroke-width="1"/>',
        f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#333" stroke-width="1"/>',
        f'<text x="{ml-12}" y="{mt-10}" text-anchor="end" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#444">{ylabel}</text>',
        f'<text x="{ml+pw/2}" y="{height-12}" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#444">{xlabel}</text>',
    ]
    for xv in _nice_ticks(x0, x1, 6):
        xx = X(xv)
        parts.append(
            f'<line x1="{xx}" y1="{mt}" x2="{xx}" y2="{mt+ph}" stroke="#eee" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{xx}" y="{mt+ph+16}" text-anchor="middle" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#555">'
            f"{xv:.3g}</text>"
        )
    for yv in _nice_ticks(y0, y1, 5):
        yy = Y(yv)
        parts.append(
            f'<line x1="{ml}" y1="{yy}" x2="{ml+pw}" y2="{yy}" stroke="#eee" stroke-width="1"/>'
        )
        lab = f"{yv:.0f}" if yv >= 10 else f"{yv:.1f}"
        parts.append(
            f'<text x="{ml-8}" y="{yy+3}" text-anchor="end" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#555">'
            f"{lab}</text>"
        )

    legend_y = mt + 14
    legend_x = ml + pw - 8
    for name, pts, color, style in series:
        if not pts:
            continue
        if style in ("line", "both"):
            d = " ".join(
                f"{'M' if i == 0 else 'L'}{X(x):.2f},{Y(y):.2f}"
                for i, (x, y) in enumerate(pts)
            )
            parts.append(
                f'<path d="{d}" fill="none" stroke="{color}" stroke-width="1.8"/>'
            )
        if style in ("markers", "both"):
            for x, y in pts:
                if "exp" in name.lower() or name == "Experiment":
                    # squares for experiment
                    s = 5.0
                    parts.append(
                        f'<rect x="{X(x)-s:.2f}" y="{Y(y)-s:.2f}" width="{2*s}" height="{2*s}" '
                        f'fill="{color}" stroke="white" stroke-width="1"/>'
                    )
                else:
                    parts.append(
                        f'<circle cx="{X(x):.2f}" cy="{Y(y):.2f}" r="4.2" '
                        f'fill="{color}" stroke="white" stroke-width="1"/>'
                    )
        legend_y += 16
        parts.append(
            f'<line x1="{legend_x-78}" y1="{legend_y}" x2="{legend_x-54}" y2="{legend_y}" '
            f'stroke="{color}" stroke-width="2"/>'
        )
        parts.append(
            f'<text x="{legend_x-48}" y="{legend_y+3}" text-anchor="start" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#222">'
            f"{name}</text>"
        )

    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_svgs(rows: list[dict]) -> tuple[Path, Path]:
    exp_Q = [(r["m_dot"], r["Q_exp_W"]) for r in rows if r["Q_exp_W"] is not None]
    cfd_Q = [(r["m_dot"], r["Q_CFD_W"]) for r in rows if r["Q_CFD_W"] is not None]
    exp_P = [(r["m_dot"], r["dP_exp_Pa"]) for r in rows if r["dP_exp_Pa"] is not None]
    cfd_P = [(r["m_dot"], r["dP_CFD_Pa"]) for r in rows if r["dP_CFD_Pa"] is not None]

    p_q = OUT_FIG / "validation_heat_transfer_rate.svg"
    p_p = OUT_FIG / "validation_pressure_drop.svg"
    svg_xy(
        [
            ("Experiment (Fig.5)", exp_Q, "#111111", "both"),
            ("OpenFOAM", cfd_Q, "#c45c26", "both"),
        ],
        p_q,
        title="Fig. 5(a) — heat transfer rate vs mass-flow (SG-STHX)",
        xlabel="ṁ [kg/s]",
        ylabel="Q [W]",
    )
    svg_xy(
        [
            ("Experiment (Fig.5)", exp_P, "#111111", "both"),
            ("OpenFOAM", cfd_P, "#2f6fed", "both"),
        ],
        p_p,
        title="Fig. 5(b) — shell-side pressure drop vs mass-flow (SG-STHX)",
        xlabel="ṁ [kg/s]",
        ylabel="Δp [Pa]",
    )
    return p_q, p_p


def main() -> None:
    exp_q, exp_dp = parse_fig5_experiment(REF_MD)
    print(f"Fig.5 experiment from {REF_MD.name if REF_MD.is_file() else 'fallback'}")
    for m in M_DOTS:
        print(f"  m={m:.3f}: Q_exp={exp_q[m]:.0f} W  dP_exp={exp_dp[m]:.0f} Pa")

    cfd = load_cfd_points()
    rows = build_table(cfd, exp_q, exp_dp)
    paper_csv = write_paper_csv(exp_q, exp_dp)
    csv_path = write_csv(rows)
    md_path = write_md(rows)
    fig_q, fig_p = write_svgs(rows)
    print(f"wrote {paper_csv.relative_to(CASE)}")
    print(f"wrote {csv_path.relative_to(CASE)}")
    print(f"wrote {md_path.relative_to(CASE)}")
    print(f"wrote {fig_q.relative_to(CASE)}")
    print(f"wrote {fig_p.relative_to(CASE)}")
    for r in rows:
        if r["Q_CFD_W"] is not None:
            print(
                f"  OF m={r['m_dot']:.3f}: Q={r['Q_CFD_W']:.1f} W "
                f"({r['err_Q_vs_exp_pct']:+.1f}% vs exp), "
                f"dP={r['dP_CFD_Pa']:.1f} Pa ({r['err_dP_vs_exp_pct']:+.1f}% vs exp)"
            )


if __name__ == "__main__":
    main()
