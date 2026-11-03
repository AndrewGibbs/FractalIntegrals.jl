# Changes to existing code

Every AI change to code that already existed is recorded here (ground rule 3). See `git log` for the exact diffs, and `dev_notes/bugs/bugs.ipynb` for what each bug did.

## Bug fixes, approved by Andrew on 2026-11-03

Tests are in `test/regression_tests.jl`.

- **Block operators.** `polygonapprox/polygonapprox.jl`: `function ∪` → `function Base.:∪`. It created a new `∪` inside the package, so internal calls failed. Broken since `4f26c84`.

- **3D Helmholtz operator.** `operators/operatorpresets.jl`: `ComplexF64{R}` → `Complex{R}`. Broken since `0e0055f`.
