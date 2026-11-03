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

- **Laplace potential.** `potentials/potentialpresets.jl`, `operators/kernels.jl`: `singlelayer_potential_laplace` did not use the operator's Green's function. It used `n = ambient_dimension + 1`, so in 2D the kernel was `1/(2πr)` and in 3D it was `1/(8πr²)`. It now uses `-log r/(2π)` and `1/(4πr)`, with screen versions for N=1 in 2D and N=2 in 3D, as the Helmholtz potential does. Checked against the exact solution on a segment, the capacitance of a square plate, and `-Δu = φ` on a square (`test/exact_solution_tests.jl`). The operator itself was correct. The `gamma` import from SpecialFunctions was only used here, so is removed.

- **Block operators.** `discretise` on a `BlockOperator` errored, because it expected one measure and one basis. It now builds one mesh and quadrature rule per component, and `\` and the potentials work on the result. Collocation gives an `ArgumentError`, as block collocation is a separate to-do item. Along the way:
  - `discretisation/sesquilinearforms.jl`: off-diagonal blocks were missing the normalisation factors. **Any matrix built by hand with `get_galerkin_matrix` on a tuple of bases changes.**
  - `operators/operatorpresets.jl`: `ambient_dimension` was not passed to the diagonal blocks, so the 3D kernel gave a 2D diagonal. `check_ambient_dimension` is now called, since `ambient_dimension = 1` errored off the diagonal. The default in `dom2codom_*` could not run.
  - `discretisation/galerkin.jl`: `reps` and other keywords no longer reach the off-diagonal blocks, which can't take them. The matrix type comes from its blocks rather than always being complex.
  - `fractalgeometry/binaryoperations.jl`: `embed_into_same_dimension` compared with `fill(x) == tuple`, which is always false, so it re-embedded when it had nothing to do.
  - `discretisation/basis.jl`, `potentials/potentialpresets.jl`: `get_ambient_dimension` for bases, so that presets work on a `UnionBasis`.
  - Check: splitting the Cantor set or dust into its pieces gives the same matrix as the single operator, to about 1e-16.

- **Docs references.** `docs/src/integralequations.md`: the `@docs` blocks now use `FractalIntegrals.singlelayer_operator_laplace` and `_helmholtz`, as these aren't exported. Without the prefix the docs build fails, and so did the Make docs CI job. This is the one edit to docs in this round.

- **`SingularIntegralOperator`.** `operators/integraloperators.jl`: its type parameter `Z` wasn't in any field, so it couldn't be constructed without naming all its parameters, and the `@hausdorffdefault` method always errored. It now has a constructor, which works under collocation. Galerkin gives an `ArgumentError`, because the singular algorithm needs the split into `r^-s` and a smooth part: use `SeparableIntegralOperator`, or collocation.

- **`I + S` under collocation.** `discretisation/collocation.jl`: added `compute_col_entry` for `IdentityOperator`, `ScaledOperator` and `SumOperator`, so a second-kind equation can be collocated. `discretisation/discretisation.jl`: added their default mesh widths, so `discretise` needs no `h_mesh`. `discretisation/galerkin.jl`: their Galerkin methods now accept keywords such as `reps`, as the integral operators do.

- **Numeric types.** Presets, `s_energy` and Laplace Galerkin now work in Float32 and BigFloat, and keep that type. Float64 results change only by rounding (`x/(2π)` became `x/2/π`).
  - `fractalgeometry/invariantmaps.jl`: `InvariantMap` always made Float64, and the rotation loop used a Float64 range. `fractalgeometry/attractors.jl`: `dimH`.
  - `fractalpresets/`: `1/3` and `T(1/3)` became `1/T(3)`, and so on. This wasn't only cosmetic: `sierpinskicarpet(Float32)` failed its own assertion, and `vicsek` was silently disconnected for anything but Float64, because `ρ == 1/3` was false.
  - `quadrature/senergy.jl`: `eps()` became `eps(T)` in the similarity tests, which made Float32 tests far too strict, and arrays are now typed. `quadrature/jacobimatrices.jl`, `discretisation/sesquilinearforms.jl`, `discretisation/basis.jl`, `discretisation/collocation.jl`, `operators/operatorpresets.jl`, `operators/kernels.jl`: Float64 constants and arrays.
  - `quadrature/quadrature.jl`: the default rule is barycentre if the type isn't a LAPACK one, as Gauss needs `eigen`. `FractalIntegrals.jl` imports `BlasFloat` for this. BigFloat Gauss would need GenericLinearAlgebra, a new dependency, so I haven't added it.
  - Not possible: BigFloat Helmholtz, because SpecialFunctions has no BigFloat `hankelh1`.

- **Symmetry centres.** `fractalpresets/compactattractors.jl`: the dihedral groups of `cantordust`, `sierpinskitriangle`, `sierpinskicarpet` and `vicsek` rotated about the origin rather than the centre of the fractal, so didn't map it to itself (points moved up to about 1.4). They now use `[1/2, 1/2]`, and for the triangle `[1/2, √3/6]/(2(1-ρ))` with an angle offset of π/6. Results on these four presets are unchanged: energies and Galerkin matrices match exactly, before and after. Andrew approved this on 2026-11-03.

- **Gauss rule on measures of mass other than one.** `quadrature/jacobimatrices.jl`: `A[1]` was multiplied by `μ.suppmeasure`, though the rest of the recurrence is for the normalised measure, so the Jacobi matrix had NaNs and `eigen` threw. `quadrature/gauss.jl`: the weights now sum to μ(Γ). Results on measures of mass one are unchanged. Sub-measures such as `μ[1]` now work with the default quadrature.

- **Vicsek connectedness.** `fractalpresets/compactattractors.jl`: the corner pieces touch the centre, but only one triangle of the matrix said so. So the pairs `([m],[5])` were treated as smooth integrals in the singular algorithm. **Numbers on the Vicsek fractal change,** by about 1e-6 in energies (s = 0.5).

- **Choice of integration method.** `discretisation/sesquilinearforms.jl`, `quadrature/senergy.jl`: Galerkin entries were sent to the singular method only if the bounding balls overlapped, `dist⁻ <= 0`. For touching pieces the balls can just touch, so rounding decided, and some touching pairs got plain quadrature across their singularity. The test is now `connected_by_matrix`, followed by the existing check for a similar canonical integral, whose list comes from the connectedness matrix. Andrew asked for connectivity to decide. **Numbers change** on fractals with touching pieces, by about 1e-5 in the matrix on Sierpinski and Vicsek. Assembly is slower where most pieces touch, since more pairs reach the similarity check: about 9 times on Sierpinski (0.06 s to 0.55 s, n = 243) and 7 times on Vicsek (n = 625), but unchanged on the dust and Koch curve. Collocation still uses the balls.
