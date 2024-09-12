```@setup tutorial
using FractalIntegrals
using Plots; gr()
Plots.reset_defaults()
```

# Constructing Fractals and Fractal measures

This package focuses on self-similar fractals, which can be described by a finite set of self-similar maps. These maps are called *similarities*.

## Similarities
Each similarity is of the form:

```math
s_m(x)=\rho_mA_mx + \delta_m,\quad x\in\mathbb{R}^n,
```

where ``\rho_m\in(0,1)`` is the contraction factor, ``A_m\in\R^{N\times N}`` is a rotation matrix, and ``\delta\in\R^{N}`` is a translation vector. Here we call $N\in\mathbb{N}$ the ambient dimension.

```@docs
Similarity
```

## Attractors

An iterated function system (IFS) is a set of $M$ similaritites ``\{s_m\}_{m=1}^M``, and an *attractor* of an iterated function system $\Gamma$ satisfies

```math
\Gamma  = \bigcup_{m=1}^M s_m(\Gamma)
```

```@docs
Attractor
```

In general, an IFS attractor has a non-integer [Hausdorff dimension](https://en.wikipedia.org/wiki/Hausdorff_dimension). The Hausdorff dimension $d$ satisfies
```math
\sum_{m=1}^M\rho_m^d = 1.
```
The following example constructs a Cantor set and displays its dimension.

```@example tutorial
s₁ = Similarity(1/3,0)
s₂ = Similarity(1/3,2/3)
Γ = Attractor(s₁, s₂)
print(Γ)
```
```@example tutorial
println("Full Hausdorff dimension of Cantor Set:", Γ.d)
```

### Plotting attractors
`Plots.plot` has a method for a `Attractor` type.
```@example tutorial
plot(Γ)
```
### Presets

You can create your own attractor from scratch, as explained above. But a range of preset attractors are available.

```@docs
getfractal
```

```@example tutorial
Γ = getfractal("sierpinskitriangle")
plot(Γ)
```

### Sub-components

It is often convenient to talk about a subcomponent of a fractal $\Gamma$ which is a scaled copy of the original $\Gamma$. This is particularly useful when meshing the fractal with fractal mesh elements. To do this, vector indexing is used, for example for $\mathbf{m}=[m_1,\ldots,m_\ell]$,

```math
\Gamma_{\mathbf{m}} = s_{m_1}\circ \ldots \circ s_{m_\ell} (\Gamma)
```

This can be achieved by treating attractors as vectors.

```@example tutorial
Γ₁ = Γ[1]
𝐦 = [3,2,3,1]
Γ₃₂₃₁ = Γ[𝐦]
plot!(Γ₁,markersize=1)
plot!(Γ₃₂₃₁,markersize=1.5)
```

## Measures

We now consider a measure $\mu$ supported on a fractal attractor $\Gamma$. This is necessary to define integrals and integral equations on $\Gamma$. For a attractor defined with $M$ similarities, an invariant measure $\mu$ is defined by an associated set of *probabilities* $p_1,\ldots,p_M$ satisfying $\sum_{m=1}^M p_m =1$, such that

```math
\mu(\Gamma_\mathbf{m}) = \left(\prod_{i=1}^\ell p_{m_i}\right)\mu(\Gamma).
```

```@docs
InvariantMeasure
```

The *Hausdorff Measure* $\mathcal{H}^d$ is a special case of an *Invariant Measure*, where $p_m=\rho_m^d$ for $m=1,\ldots,M$. This ensures that the mass of the measure is, in some sense, spread evenly across $\Gamma$. A useful consequence is that if $\Gamma$ is invariant under any non-trivial rotations/reflections, then this is inherited by the Hausdorff measure.

```@docs
HausdorffMeasure
```

The Hausdorff measure is the most natural measure for any fractal, and is used by default throughout `FractalIntegrals`. Recall that from an above example, $\Gamma$ represents the Sierpinski triangle.

```@example tutorial
# no weights specified
𝓗ᵈ = InvariantMeasure(Γ)
print(typeof(𝓗ᵈ))
```
### Plotting measures
When `plot` is used on an `InvariantMeasure`, colouring is used to represent the distribution of mass.

```@example tutorial
# create random set of probabilites
𝐩 = rand(3)
𝐩 = 𝐩/sum(𝐩)
μ = InvariantMeasure(Γ, 𝐩)
plot(μ, markersize = 1)
```