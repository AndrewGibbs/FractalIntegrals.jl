# Changes to existing code

Every AI change to code that already existed is recorded here (ground rule 3). See `git log` for the exact diffs, and `dev_notes/bugs/bugs.ipynb` for what each bug did.

## Bug fixes, approved by Andrew on 2026-11-03

Tests are in `test/regression_tests.jl`.

- **Block operators.** `polygonapprox/polygonapprox.jl`: `function ∪` → `function Base.:∪`. It created a new `∪` inside the package, so internal calls failed. Broken since `4f26c84`.

- **3D Helmholtz operator.** `operators/operatorpresets.jl`: `ComplexF64{R}` → `Complex{R}`. Broken since `0e0055f`.

- **Preset dimensions.** `fractalpresets/compactattractors.jl`: Sierpinski `log(5)` → `log(3)`; carpet `log(3)` → `log(1/ρ)`. **This changes numerical results on the Sierpinski triangle**, and on carpets with ρ ≠ 1/3.

- **Barycentre weights.** `quadrature/barycentrerule.jl`: the attractor type was in the position of the parameter `N`, so the method never applied.

- **`getmeasure`.** `fractalpresets/measurepresets.jl`: after its warning, set `suppmeasure = one(T)`, as the warning describes. It stayed as the type `Nothing`.
