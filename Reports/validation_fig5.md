# Model validation — Zahid Fig. 5 (SG-STHX)

Primary reference: **experimental** points digitised in [MAIN REFERENCE](../../../Reference/(MAIN%20REFERENCE)%20Experimental%20and%20CFD%20Simulation%20Study%20of%20Shell%20and%20Tube%20Heat%20Exchangers%20with%20Different%20Baffle%20Segment%20Configurations.md) Fig. 5.

**OpenFOAM:** `m037` counterflow, k–ε, ṁ = 0.037 kg/s, Tin = 371.47 / 298 K (Fig. 6 legend max/min as inlets).

Also: our **Kern** (Q) / **Esso** (Δp) closed-form rebuild with the same geometry assumptions (pitch 12 mm, Tin from Fig. 6). Paper Fig. 5 Kern/Esso chart curves were not digitised.

## Validation vs experiment

| ṁ [kg/s] | Q<sub>exp</sub> [W] | Q<sub>OF</sub> [W] | err Q [%] | Δp<sub>exp</sub> [Pa] | Δp<sub>OF</sub> [Pa] | err Δp [%] |
|---:|---:|---:|---:|---:|---:|---:|
| 0.027 | 2230 | — | — | 2000 | — | — |
| 0.032 | 2300 | — | — | 2900 | — | — |
| 0.037 | 2400 | 2304.7 | -4.0 | 3750 | 3396.8 | -9.4 |
| 0.042 | 2520 | — | — | 4750 | — | — |
| 0.047 | 2650 | — | — | 5800 | — | — |

## Validation vs Kern / Esso (same assumed Tin & pitch)

Paper claims CFD within ~15% of Kern (Q) and ~10% of Esso (Δp). Errors below are **OpenFOAM − closed-form** / closed-form.

| ṁ [kg/s] | Q<sub>Kern</sub> [W] | Q<sub>OF</sub> [W] | err Q [%] | Δp<sub>Esso</sub> [Pa] | Δp<sub>OF</sub> [Pa] | err Δp [%] |
|---:|---:|---:|---:|---:|---:|---:|
| 0.027 | 1021.5 | — | — | 21.6 | — | — |
| 0.032 | 1138.2 | — | — | 29.3 | — | — |
| 0.037 | 1247.3 | 2304.7 | +84.8 | 38.2 | 3396.8 | +8801.0 |
| 0.042 | 1350.1 | — | — | 48.0 | — | — |
| 0.047 | 1447.5 | — | — | 58.8 | — | — |

## Figures

- [`figures/validation_heat_transfer_rate.svg`](figures/validation_heat_transfer_rate.svg)
- [`figures/validation_pressure_drop.svg`](figures/validation_pressure_drop.svg)

## Notes

- Errors are **OpenFOAM − reference** / reference.
- Tin from Fig. 6 colour-bar extremes: **371.47 K** (hot), **298 K** (cold).
- Shell nozzles: equivalent Ø **5.5 mm** (D_shell/8), uniform W×L stub into the shell (multi-cell in Y).
- Our Kern/Esso rebuild uses the same Tin/pitch assumptions — relative check only; not the paper’s published Kern/Esso curves.
