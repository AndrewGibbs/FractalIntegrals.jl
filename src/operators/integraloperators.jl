abstract type IntegralOperator{M<:Measure, Z<:Number} <: FractalOperator{M} end
abstract type AbstractSingularIntegralOperator{M <: AbstractInvariantMeasure, Z} <: IntegralOperator{M, Z} end
abstract type AbstractSeparableIntegralOperator{M, Z} <: AbstractSingularIntegralOperator{M, Z} end

"""
    SeparableIntegralOperator{M, F1, F2, S, Z} <: AbstractSeparableIntegralOperator{M, Z}

An integral operator whose kernel splits into a singular part and a smooth (Lipschitz)
part, of the form

```math
K(x,y) = K_{\\mathrm{smooth}}(x,y) + c \\cdot K_s(x,y),
```

where ``K_{\\mathrm{smooth}} \\in C^{2}(\\mathbb{R}^n)`` and

```math
K_s(x,y) := \\begin{cases} \\log|x-y|, & s=0, \\\\ |x-y|^s, & s>0, \\end{cases}
```

for some ``s \\geq 0``. This decomposition enables efficient numerical treatment:
the smooth part is integrated by standard fractal quadrature rules, and the singular
part is handled by a dedicated singular integration algorithm.

This structure is appropriate for Laplace-type kernels (no oscillation). For
oscillatory kernels (e.g. Helmholtz), see `OscillatorySeparableIntegralOperator`.
For a purely singular kernel with no smooth part, see `SingularIntegralOperator`.

# Fields
- `measure :: M`: The invariant measure ``\\mu`` supported on ``\\Gamma``.
- `kernel :: F1`: The full kernel function ``K(x,y)``.
- `lipschitzpart :: F2`: The smooth part ``K_{\\mathrm{smooth}}(x,y)`` of the kernel.
- `s :: S`: The singularity exponent (``s=0`` for logarithmic, ``s>0`` for algebraic).
- `singularconst :: Z`: The constant ``c`` scaling the singular part.
- `symmetric :: Bool`: Whether the operator is self-adjoint.
"""
struct SeparableIntegralOperator{
        M <: AbstractInvariantMeasure,
        F1 <: Function,
        F2 <: Function,
        S <: Real,
        Z <: Number
        } <: AbstractSeparableIntegralOperator{M, Z}
    measure::M
    kernel::F1
    lipschitzpart::F2
    s::S
    singularconst::Z
    symmetric::Bool   
end

struct SmoothIntegralOperator{
        Mdom <: AbstractInvariantMeasure,
        Mcodom <: AbstractInvariantMeasure,
        F <: Function,
        Z <: Number
        } <: IntegralOperator{Mdom, Z}
    measure :: Mdom
    measure_codom :: Mcodom
    kernel :: F
    symmetric :: Bool  
end

"""
    SingularIntegralOperator{M, F, Z} <: AbstractSingularIntegralOperator{M, Z}

An integral operator with a purely singular kernel, i.e. with no additional smooth part.
The kernel takes the form

```math
\\mathcal{K}\\,\\phi\\,(x) = \\int_\\Gamma K(x,y)\\,\\phi(y)\\,\\mathrm{d}\\mu(y),
\\quad x \\in \\Gamma,
```

where ``K`` is singular (e.g. a Laplace-type energy kernel). Unlike
`SeparableIntegralOperator`, no kernel splitting is assumed; the entire kernel is
treated as singular.

# Fields
- `measure :: M`: The invariant measure ``\\mu`` supported on ``\\Gamma``.
- `kernel :: F`: The singular kernel function ``K(x,y)``.
- `symmetric :: Bool`: Whether the operator is self-adjoint.
"""
struct SingularIntegralOperator{
            M <: AbstractInvariantMeasure,
            F <: Function,
            Z <: Number
        } <: AbstractSingularIntegralOperator{M, Z}
    measure :: M
    kernel :: F
    symmetric :: Bool  
end

# ultimately, should bin off the below type, its too messy
struct OscillatorySeparableIntegralOperator{
        M <: AbstractInvariantMeasure,
        F1 <: Function,
        F2 <: Function,
        S <: Real,
        Z <: Number,
        } <: AbstractSeparableIntegralOperator{M, Z}
    measure::M
    kernel::F1
    lipschitzpart::F2
    s::S
    singularconst::Z
    symmetric::Bool
    wavenumber::S
end

@hausdorffdefault SingularIntegralOperator
@hausdorffdefault OscillatorySeparableIntegralOperator