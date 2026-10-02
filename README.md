# Zahid lab SG-STHX — ṁ = 0.037 kg/s

Variation of the Table-1 CHT case. Both tube and shell inlets use
`flowRateInletVelocity` with `massFlowRate 0.037;`.

```sh
./Allclean && ./buildMesh   # if mesh/BCs need refresh
./Run
```

Parent runner: [`../Run-variations`](../Run-variations).
