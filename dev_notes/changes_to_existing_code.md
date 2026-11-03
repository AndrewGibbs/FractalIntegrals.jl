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

## Formatting and merge, approved by Andrew on 2026-11-03

- **All Julia files reformatted.** Andrew asked to merge Zoïs Moitier's work from `origin/main`, which ran JuliaFormatter over the package. It also rewrites code forms: adds explicit `return`, turns one-line functions into `function ... end`, and `?:` into `if/else`. Our side was formatted the same way first. The merged code parses identically to our formatted code, and all tests pass.
- **Comments kept.** Zoïs's formatting run removed some comments, such as `# self-adjoint`. Where the two sides conflicted the merge kept our version, so these remain.
- **Also from that merge:** CI tests on Julia 1.10, latest and nightly (no longer 1.7.2), `test/Manifest.toml` and a notebook checkpoint removed, and a larger `.gitignore`.

## Julia version, approved by Andrew on 2026-11-03

- **`Project.toml`:** `julia = "1.7 - 1.12"` → `julia = "1.10"`, meaning 1.10 or later. The upper cap appeared to stop the package installing on the latest Julia and on nightly in CI, and 1.10 is now the oldest version CI tests.

## Second round of fixes, approved by Andrew on 2026-11-03

- **Laplace potential.** `potentials/potentialpresets.jl`, `operators/kernels.jl`: `singlelayer_potential_laplace` did not use the operator's Green's function. In 2D it evaluated `-2π/0 × log r`-style nonsense through `1/((n-2)Aₙ)`, and in 3D the constant was wrong and it used `Aₙ` for the wrong dimension. It now uses `-log r/(2π)` and `1/(4πr)`, with screen versions for N=1 in 2D and N=2 in 3D, as the Helmholtz potential does. Checked against the exact solution on a segment, the capacitance of a square plate, and `-Δu = φ` on a square (`test/exact_solution_tests.jl`). The operator itself was correct. The `gamma` import from SpecialFunctions was only used here, so is removed.
