```@setup tutorial
using FractalIntegrals
using Plots; gr()
Plots.reset_defaults()
```

# Integrals on fractals

In [the previous section](makeIFS.md) we saw how to define a fractal measure $\mu$ supported on some fractal attractor $\Gamma$. Equipped with this measure, we can define the integral

```math
I_\mu[f] := \int_\Gamma~ f(x)~ \mathrm{d}\mu(x),
```
for some function $f$.

## Quadrature rules

In general, $I$ cannot be computed exactly. We must use *quadrature rules* to approximate $I$; these consist of nodes $x_j$ and weights $w_j$ for $j=1,\ldots,N$ such that
```math
I_\mu[f] \approx Q_\mu[f]:=\sum_{j=1}^N f(x_j)~ w_j
```

for non-fractal domains and measures which are sufficiently simple, a range of quadrature rules are available. See for example, [`FastGaussQuadrature`](https://juliaapproximation.github.io/FastGaussQuadrature.jl/stable/) or [`QuadGK`](https://juliamath.github.io/QuadGK.jl/stable/).

In recent decades, quadrature rules have been developed for fractal measures $\mu$. These have been implemented in `FractalIntegrals`.

>[!NOTE]
>Contrary to classical quadrature rules, the nodes typically lie outside of $\Gamma$. Therefore, it is necessary to assume that $f$ is supported on $\mathrm{Hull}(\Gamma)$.

### Gauss Quadrature
On standard domains, Gaussian qudarature rules are the most widely used for smooth $f$. The classical process, which involves constructing orthogaonl polynomials and a Jacobi matrix, was generalised in [Ma:96](@cite) to fractal measures $\mu$ where 
```math
\mathrm{supp}(\mu) = \Gamma \subset \mathbb{R}.
```

```@docs
FractalIntegrals.gauss_quadrule
```

### Barycentre rule

Despite its simplicity, the Barycentre rule was introduced relatively recently in [GiHeMo:23](@cite). It can be interpreted as a generalisation of the midpoint rule to non-standard domains.

The idea behind this quadrature rule is to subdivide the fractal $\Gamma = \mathrm{supp}(\mu)$ into self-similar subcomponents no wider than some parameter $h>0$. A node is allocated for each subcomponent $\Gamma_{\mathbf{m}}$; this node is chosen to be the barycentre

```math
\int_{\Gamma_{\mathbf{m}}} x ~\mathrm{d}\mu(x)~/~\mu(\Gamma),
```

which ensures an $O(h^2)$ convergence rate. The reason for this is analagous to the reason that the midpoint rule converges at $O(h^2)$, because the linear errors above and below the midpoint cancel. 

```@docs
FractalIntegrals.barycentre_quadrule
```

```@example tutorial
# define invariant measure on Koch snowflake
Γ = getfractal("koch snowflake")

# get some random probability weights
p = abs.(0.5.+randn(7)/5)
# rescale the first weight due to larger size of central component
p[1] = p[1]*sqrt(3)
p = p/sum(p)

# define the invariant measure
μ = InvariantMeasure(Γ, p)

# get barycentre rule quadrature weights and nodes
h_bary = 0.01
x_bary, w_bary = FractalIntegrals.barycentre_quadrule(μ, h_bary)

# plot the distribution of weights over the snowflake
plot([x_[1] for x_ in x_bary],[x_[2] for x_ in x_bary],
    linewidth=0,
    markersize=1,
    markershape=:circle,
    markerstrokewidth=0,
    marker_z = log.(w_bary),
    axis=false,
    grid=false,
    legend = false,
    aspect_ratio=1,
    ylim=(-0.6,0.6),
    xlim=(-0.875,0.875),
    size = (875, 875),
    background=nothing)
```

### Chaos game quadrature

*Chaos game quadrature* [FoMeVr:98](@cite) is a non-deterministic approach, analagous to Monte Carlo rules. The algorithm takes an initial guess $x_0\in\mathbb{R}^n$ and repeatedly applies similarities $s_m$ at random to construct further points. Mathematically, this is:

```math
\mathbb{P}(x_{j} = s_m(x_{j-1})) = p_m,\quad\text{for }m=1,\ldots,M.
```

The random process above is repeated for $j=1,\ldots,N$; the weights are simply $1/N$.

```@docs
FractalIntegrals.chaos_quadrule
```

```@example tutorial

# get the chaos game points
n_chaos = length(x_bary)
x_chaos, w_chaos = FractalIntegrals.chaos_quadrule(μ, n_chaos)

# plot the distribution of points
plot([x_[1] for x_ in x_chaos],[x_[2] for x_ in x_chaos],
    linewidth=0,
    markersize=0.7,
    markershape=:circle,
    markerstrokewidth=0,
    axis=false,
    grid=false,
    legend = false,
    aspect_ratio=1,
    ylim=(-0.6,0.6),
    xlim=(-0.875,0.875),
    size = (875, 875),
    background=nothing)
```

## Tensor product quadrature

## Singular integrals

The quadrature rules described [above](#quadrature-rules) require some smoothness of the integrand $f$ to be accurate. We now consider a class of singular integrals:

```math
\mathcal{E}_s[\mu] := \left\{\begin{array}{cc}
    \displaystyle\int_\Gamma\int_\Gamma~ |x-y|^{-s}~ \mathrm{d}\mu(y)\mathrm{d}\mu(x),& s\neq0\\&\\
    \displaystyle\int_\Gamma\int_\Gamma~ \log|x-y|~ \mathrm{d}\mu(y)\mathrm{d}\mu(x),& s=0
    \end{array}
    \right.
```

These integrals are sometimes referred to as *s-energy* [Fa:14; Section 4.3](@cite) or the *generalised electrostatic potential* [MaVa:07; Definition 4](@cite). They are of significant interest in their own right, but also arise in integral equations posed on fractals, thus forming a key part of our integral equation solver.

```@docs
FractalIntegrals.s_energy
```

The algorithm, which was introduced in [GiHeMa:24](@cite), exploits self-similarity of $\Gamma$ and scaling properties of the singular integrand. It also exploits any rotational and reflective symmetry of the measure $\mu$, which typically only occurs when $\mu$ is the Hausdorff measure.