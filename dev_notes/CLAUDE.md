# FractalIntegrals.jl

Julia package for evaluating integrals and solving integral equations on self-similar fractals (IFS attractors). It is a rewrite of [IFSintegrals](https://github.com/AndrewGibbs/IFSintegrals). Developed by Andrew Gibbs and Josh; plans and context are in `dev_notes/next_steps.md`, and the to-do list is in `dev_notes/todo.md`.

## Ground rules

@GROUND_RULES.md

## What the code does

**Workflow:** define a fractal → put a measure on it → build an operator → discretise → solve → evaluate fields.

```julia
Γ   = getfractal("cantordust")                                 # attractor
Sₖ  = FractalIntegrals.singlelayer_operator_helmholtz(Γ, k)    # defaults to Hausdorff measure
Sₖₕ = FractalIntegrals.discretise(Sₖ, h_mesh = 0.05, h_quad = 0.01)
ϕ   = Sₖₕ \ f                                                  # Projection (basis + coefficients)
slp = FractalIntegrals.singlelayer_potential_helmholtz(ϕ, k)   # near field
ffp = FractalIntegrals.farfield_pattern(ϕ, k)                  # far field
```

Most functions are not exported, so call them as `FractalIntegrals.name`.

### Source layout (`src/`, in include order)

| Folder | Contents |
|---|---|
| `fractalgeometry/` | `Similarity` (s(x) = δ + ρAx) and `OneDimensionalSimilarity`; attractor types; measures; bounding balls and diameter; dihedral/trivial symmetry groups; translation and embedding; attractor and measure unions |
| `vectorindices/` | `VectorIndex{M,T}`: compact integer encoding of 𝐦 = [m₁,…,m_ℓ]; `subdivide_indices` (mesh) and `grade_mesh` (graded mesh) |
| `quadrature/` | Barycentre rule, chaos game, Gauss (Jacobi matrix, Mantica 1996, N=1 only), tensor products, and the s-energy singular integral algorithm |
| `operators/` | Kernels and `vecdist`; integral operator types; `singlelayer_operator_{laplace,helmholtz}`; `IdentityOperator`, `ScaledOperator`, `SumOperator`, `BlockOperator` |
| `discretisation/` | P0 basis on a quasi-uniform mesh; Galerkin and collocation; repeated-entry maps; `discretise`, `Projection`, and `\` |
| `potentials/` | `Potential`: single-layer near fields and far-field patterns, serial evaluation |
| `fractalpresets/` | `getfractal` (Cantor set/dust, Sierpinski triangle/carpet, Vicsek, Koch curve/snowflake, Heighway dragon) and `getmeasure` |
| `polygonapprox/` | `polygon_approx`: polygon outline of 2D attractors using Clipper |
| `plots/` | `Plots.plot` methods for attractors, measures and densities on Lebesgue-measure fractals |

### Type hierarchy

```
Fractal
├── AbstractAttractor{N,M,T}          N = ambient dim, M = number of maps, T = numeric type
│   ├── Attractor, OneDimensionalAttractor, OpenAttractor
│   └── AbstractHomogenousAttractor   (all ρₘ equal; has a ρ field)
│       └── HomogenousAttractor, OneDimensionalHomogenousAttractor, OpenHomogenousAttractor
└── AttractorUnion

Measure
├── AbstractInvariantMeasure{N,M,T,A}  A = attractor type, the 4th parameter
│   ├── InvariantMeasure               user-chosen probabilities
│   └── AbstractHausdorffMeasure       pₘ = ρₘᵈ
│       └── HausdorffMeasure, LebesgueMeasure
└── MeasureUnion

FractalOperator{M}
├── IntegralOperator{M,Z}
│   ├── AbstractSingularIntegralOperator
│   │   ├── SingularIntegralOperator
│   │   └── AbstractSeparableIntegralOperator: SeparableIntegralOperator, OscillatorySeparableIntegralOperator
│   └── SmoothIntegralOperator         (dom → codom, off-diagonal blocks)
└── IdentityOperator, ScaledOperator, SumOperator, BlockOperator
```

`Attractor(...)` is a type-unstable constructor that chooses the concrete type: `N == 1`, or `N == d` with integer `d` (open), otherwise general; each is then homogeneous or not.

### How the Galerkin discretisation works

1. `discretise` chooses a default mesh width: 20 DOFs per wavelength for oscillatory operators, `diam/5` otherwise. It also builds a quadrature rule on the parent measure. Barycentre is used if `h_quad > 0`; otherwise Gauss for N=1, barycentre for N>1. `h_quad` is rescaled because the rule is built on the parent and then mapped.
2. `construct_quasiuniform_p0basis` subdivides until each element has diam < `h_mesh`. It maps the quadrature rule onto each element (`PreQuadP0BasisElement`). Elements are L²-normalised.
3. `getsingularinfo` calls `construct_singularity_matrix`. This uses self-similarity and symmetry to reduce all singular self-interaction integrals to a finite linear system `A x = B r + L`, where `r` holds smooth integrals found by quadrature and `L` holds log-kernel corrections (Gibbs, Hewett, Moiola). `s_energy` uses the same algorithm.
4. `sesquilinearform` computes each matrix entry. An entry is singular if and only if it is similar to one of the fundamental singular integrals, in which case it returns scaled singular part + quadrature of the Lipschitz part. Otherwise it uses tensor quadrature of the full kernel. `connected_by_matrix` is only a filter that skips the similarity check for pairs the matrix shows are disjoint. The connectedness matrix tells Algorithm 1 (Gibbs, Hewett, Major) which integrals are singular, and only has to go deep enough for the algorithm to terminate, as the rest follow by self-similarity. It must be right: a wrong one gives silently inaccurate entries.
5. `get_galerkinreps` marks repeated entries (symmetry, and fractal block structure), and entries are computed with `@spawn`.

**Collocation** (`method = :collocation`): nodes are barycentres. `h_col` oversamples, in which case the solve uses truncated SVD. Near-singular entries use graded quadrature, which drops elements close to the node. `IntegralOperator`, `IdentityOperator` and their scaled and summed combinations are supported, with no block logic and no repeated-entry reuse.

### Conventions and pitfalls

- `[0]` is the index of the whole attractor (Γ₀ := Γ). Γ_𝐦 = s_{m₁} ∘ … ∘ s_{m_ℓ}(Γ).
- `VectorIndex` stores m₁ as the least significant digit. It overflows at around 60 levels.
- `suppmeasure` is μ(Γ) and defaults to 1. That is harmless for first-kind equations but **must be correct for second-kind ones** (see `docs/src/integralequations.md`).
- `@hausdorffdefault f` adds a method taking an attractor, which wraps it in `HausdorffMeasure`. `@tensorquad f` creates `tensor_f`.
- Kernels are vectorised: `kernel(x, y)` takes arrays of points. `energykernel(s, r)` is `r^-s`, or `log(r)` for `s = 0`.
- The field `diam` is an exact preset value or, when computed, the diameter of a bounding ball (an upper bound).
- Measures have the parameter order `{N, M, T, A}`. Writing `HausdorffMeasure{<:SomeAttractor}` constrains **N**, not `A`, and silently matches nothing. This mistake has already caused two bugs, including the block reuse issue below.

## Commands

```bash
julia --project=. -t auto -e 'using Pkg; Pkg.test()'     # test suite (test/runtests.jl)
julia --project=docs/ docs/make.jl                       # docs; needs docs env set up, see .github/workflows/docs.yml
```

Regression tests for fixed bugs are in `test/regression_tests.jl`. Also in `test/`: `exact_solution_tests.jl` (answers known independently of the package), `geometry_tests.jl`, `quadrature_tests.jl`, `operator_tests.jl`. A `@test_broken` marks a known unfixed bug. Plots, polygon approximation and graded meshes are still lightly tested, so a passing suite does not show an area works. Local Julia is 1.12.1. CI runs 1.10, latest stable and nightly on Linux, macOS and Windows. Code is formatted with JuliaFormatter.

## Files only Andrew has

These are gitignored, so they are not in the GitHub repo: `todo.odt` (the original of `dev_notes/todo.md`), `debug/`, `private_testing/`, `sols/`, the root `Manifest.toml`, and `.vscode/`.

## Known issues (found 2026-11-03)

The simple bugs are fixed: see `dev_notes/changes_to_existing_code.md`, and `dev_notes/bugs/bugs.ipynb` for what the first round did. These remain:

- **Fractal block reuse** (`get_galerkinreps` in `discretisation/matrixreps.jl`): a wrong type parameter means `fractaldiagblocks!` never runs. Correcting only the type is **wrong**: it gives an 8% matrix error on the Koch curve, which has rotating maps. It is correct to rounding on the non-rotating Cantor set and dust. On Sierpinski and Vicsek the difference was about 1e-6, measured before the bounding ball test was removed, so it needs measuring again. Proposed for `dev-speedup`, non-rotating attractors only.
- **Bounding balls** (`fractalgeometry/boundingball.jl`): the ball of radius `diam/2` does not always contain the attractor (Sierpinski reaches 1.15× the radius, Koch curve about 1.02× at its endpoints). Galerkin no longer relies on them, as self-similarity chooses the method. Collocation still does, for points near an element: measured on Sierpinski, at most about 0.07% in solutions. This needs solving before H-matrices rely on these balls.
- **Assembly time:** removing the bounding ball test made assembly about 7 to 9 times slower where most pieces touch, at coarse quadrature (see the changes file). A cheap test ahead of the similarity check is a candidate for `dev-speedup`.
- **Log capacity near one:** on the logarithmic kernel in 2D, the single layer operator is singular when the set's log capacity is one, so some separations of a union give huge densities (the Cantor set with the dust 2 above it has condition 564 and total mass about -600). This is the mathematics rather than a bug, but there is no warning.

**Import time** (`dev_notes/next_steps.md` asks why it is slow): with a warm cache it is fast (`using Plots` 0.9 s, then `using FractalIntegrals` 0.2 s). A fresh `Pkg.test` environment precompiled 201 dependencies in 134 s, including about 59 s for Plots and 3 s for FractalIntegrals. So precompiling dependencies, mostly Plots, is the likely cause of slow imports, but what triggers it on Andrew's machine has not been investigated.
