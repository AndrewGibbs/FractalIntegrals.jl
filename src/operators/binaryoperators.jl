"""
    IdentityOperator{M} <: FractalOperator{M}

Represents the identity operator ``\\mathcal{I}`` in an integral equation on a fractal
measure support ``\\Gamma``, defined by

```math
\\mathcal{I}\\,\\phi\\,(x) = \\phi(x), \\quad x \\in \\Gamma.
```

This is typically used as part of a second-kind integral equation of the form

```math
(\\lambda\\,\\mathcal{I} + \\mathcal{K})\\,v = f, \\quad \\text{on } \\Gamma,
```

where ``\\lambda \\in \\mathbb{C}`` and ``\\mathcal{K}`` is an `IntegralOperator`.

# Fields
- `measure :: M`: The invariant measure ``\\mu`` supported on ``\\Gamma``.

# Examples
```julia
Γ = getfractal("cantorset")
μ = HausdorffMeasure(Γ)
I = IdentityOperator(μ)
```
"""
struct IdentityOperator{M<:AbstractInvariantMeasure} <:FractalOperator{M}
    measure :: M
end

@hausdorffdefault IdentityOperator

# ----------------------------------------------------------------------#

struct ScaledOperator{M<:AbstractInvariantMeasure,
                    K<:FractalOperator,
                    T<:Number
                    } <: FractalOperator{M}
    measure :: M
    operator :: K
    λ :: T
end

Base.:*(λ::Number, op::FractalOperator) = ScaledOperator(op.measure, op, λ)

Base.:*(λ::Number, op::ScaledOperator) = ScaledOperator(op.measure, op.operator, op.λ*λ)

Base.:*(op::FractalOperator, λ::Number) = λ*op


# ----------------------------------------------------------------------#

struct SumOperator{M<:AbstractInvariantMeasure,
                K<:FractalOperator,
                T<:FractalOperator
                } <: FractalOperator{M}
    measure :: M
    operator1 :: K
    operator2 :: T
end

function Base.:+(op1::FractalOperator, op2::FractalOperator)
    @assert op1.measure == op2.measure "measures of summed operators must match"
    SumOperator(op1.measure, op1, op2)
end

# ------------------------------------------------------------------------ #

struct BlockOperator{M<:MeasureUnion, B<:Matrix} <: FractalOperator{M}
    measures :: M
    operators :: B
    symmetric :: Bool
end

Base.getindex(K::BlockOperator, args...) = getindex(K.operators, args...)

#                      ,,,,,,,,,,,,,,,,,,,,,,,,,,,,............................................;[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]])