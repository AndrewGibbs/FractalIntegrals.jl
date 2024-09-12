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

```math
s_m(x)=\rho_mA_mx + \delta_m,\quad x\in\mathbb{R}^n,
```

```@docs
Attractor
```
### Plotting attractors
In the following example, we build and plot a Cantor set from scratch.
```@example tutorial
using FractalIntegrals
s₁ = Similarity(1/3,0)
s₂ = Similarity(1/3,2/3)
Γ = Attractor(s₁, s₂)
plot(Γ)
```
Note that `Plots.plot` has a method for a `Attractor` type.

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

## Fractal Measures

We now consider a measure $\mu$ supported on a fractal attractor $\Gamma$. This is necessary to define integrals and integral equations on $\Gamma$.