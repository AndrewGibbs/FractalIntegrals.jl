# Changes to existing code

Every AI change to code that already existed is recorded here (ground rule 3). See `git log` for the exact diffs, and `dev_notes/bugs/bugs.ipynb` for what each bug did.

## Bug fixes, approved by Andrew on 2026-11-03

Tests are in `test/regression_tests.jl`.

- **Block operators.** `polygonapprox/polygonapprox.jl`: `function ∪` → `function Base.:∪`. It created a new `∪` inside the package, so internal calls failed. Broken since `4f26c84`.

- **3D Helmholtz operator.** `operators/operatorpresets.jl`: `ComplexF64{R}` → `Complex{R}`. Broken since `0e0055f`.

- **Preset dimensions.** `fractalpresets/compactattractors.jl`: Sierpinski `log(5)` → `log(3)`; carpet `log(3)` → `log(1/ρ)`. **This changes numerical results on the Sierpinski triangle**, and on carpets with ρ ≠ 1/3.

- **Barycentre weights.** `quadrature/barycentrerule.jl`: the attractor type was in the position of the parameter `N`, so the method never applied.

- **`getmeasure`.** `fractalpresets/measurepresets.jl`: after its warning, set `suppmeasure = one(T)`, as the warning describes. It stayed as the type `Nothing`.

- **Single-element mesh.** `quadrature/quadrature.jl`: `mapquadrule` returns early for index `[0]`, which reached `ifs[0]`. This was the bug in the to-do list.

- **Nested inner product.** `discretisation/sesquilinearforms.jl`: compared an index to a basis element, so was always false. It also used the larger element's measure, but the overlap is the smaller element.

- **Small definition errors.** `fractalgeometry/unions.jl`: `.attractors` → `.measures`. `fractalgeometry/similarities.jl`: `sim_map_inv` used a field that doesn't exist, and divided on the wrong side. `FractalIntegrals.jl`: removed the export of `FractaPresets`, which doesn't exist.

- **`getindex`.** `operators/binaryoperators.jl`: the argument was *named* `BlockOperator`, so the method applied to every type. It now applies to `BlockOperator` only.

- **Docs and examples.** `README.md`: `Pkg.install` → `Pkg.add`. `docs/src/integralequations.md`, `operators/integraloperators.jl` and `quadrature/senergy.jl`: kernel sign, `|x-y|^s` → `|x-y|^{-s}`, to match the code. `fractalgeometry/attractors.jl`: docstring example `Bool(ones(3,3))` → `ones(Bool, 3, 3)`. `examples/two_dimensional_scattering.jl` and `.ipynb`: replaced `Γ.n`, and in the `.jl`, `ϕ`, neither of which exist there.
