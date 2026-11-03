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
4. `sesquilinearform` computes each matrix entry. If `dist⁻` of the two supports is ≤ 0, it looks for a similar canonical singular integral and returns scaled singular part + quadrature of the Lipschitz part. Otherwise it uses tensor quadrature of the full kernel.
5. `get_galerkinreps` marks repeated entries (symmetry, and fractal block structure), and entries are computed with `@spawn`.

**Collocation** (`method = :collocation`): nodes are barycentres. `h_col` oversamples, in which case the solve uses truncated SVD. Near-singular entries use graded quadrature, which drops elements close to the node. Only `IntegralOperator` is supported, with no block logic and no repeated-entry reuse.

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

Regression tests for fixed bugs are in `test/regression_tests.jl`. Many code paths are still untested, so a passing suite does not show an area works. Local Julia is 1.12.1. `Project.toml` allows 1.7–1.12 and CI runs 1.7.2, latest stable and nightly on Linux, macOS and Windows.

## Files only Andrew has

These are gitignored, so they are not in the GitHub repo: `todo.odt` (the original of `dev_notes/todo.md`), `debug/`, `private_testing/`, `sols/`, the root `Manifest.toml`, and `.vscode/`.

## Known issues (found 2026-11-03)

The simple bugs are fixed: see `dev_notes/changes_to_existing_code.md`, and `dev_notes/bugs/bugs.ipynb` for what they did. These two remain:

- **Fractal block reuse** (`discretisation/matrixreps.jl:108`): a wrong type parameter means `fractaldiagblocks!` never runs. Correcting only the type is **wrong**: it gives an 8% matrix error on the Koch curve, which has rotating maps. It is correct to rounding on the non-rotating Cantor set and dust. On Sierpinski and Vicsek the difference is about 1e-6, in line with the rounding issue below. Proposed for `dev-speedup`, non-rotating attractors only.
- **Bounding balls** (`fractalgeometry/boundingball.jl`): the ball of radius `diam/2` does not always contain the attractor (Sierpinski reaches 1.15× the radius, Koch curve about 1.02× at its endpoints), and for touching pieces `dist⁻ <= 0` is decided by rounding. Measured on Sierpinski, the effect is small: about 1e-7 in Galerkin far fields, and at most about 0.07% in collocation solutions. Left unchanged for now, but needs solving before H-matrices rely on these balls.

**Import time** (`dev_notes/next_steps.md` asks why it is slow): with a warm cache it is fast (`using Plots` 0.9 s, then `using FractalIntegrals` 0.2 s). A fresh `Pkg.test` environment precompiled 201 dependencies in 134 s, including about 59 s for Plots and 3 s for FractalIntegrals. So precompiling dependencies, mostly Plots, is the likely cause of slow imports, but what triggers it on Andrew's machine has not been investigated.
