struct IdentityOperator{M<:AbstractInvariantMeasure} <: FractalOperator{M}
    measure::M
end

@hausdorffdefault IdentityOperator

# ----------------------------------------------------------------------#

struct ScaledOperator{M<:AbstractInvariantMeasure,K<:FractalOperator,T<:Number} <:
       FractalOperator{M}
    measure::M
    operator::K
    λ::T
end

Base.:*(λ::Number, op::FractalOperator) = ScaledOperator(op.measure, op, λ)

Base.:*(λ::Number, op::ScaledOperator) = ScaledOperator(op.measure, op.operator, op.λ * λ)

Base.:*(op::FractalOperator, λ::Number) = λ * op

# ----------------------------------------------------------------------#

struct SumOperator{M<:AbstractInvariantMeasure,K<:FractalOperator,T<:FractalOperator} <:
       FractalOperator{M}
    measure::M
    operator1::K
    operator2::T
end

function Base.:+(op1::FractalOperator, op2::FractalOperator)
    @assert op1.measure == op2.measure "measures of summed operators must match"
    return SumOperator(op1.measure, op1, op2)
end

# ------------------------------------------------------------------------ #

struct BlockOperator{M<:MeasureUnion,B<:Matrix} <: FractalOperator{M}
    measures::M
    operators::B
    symmetric::Bool
end

Base.getindex(BlockOperator, args...) = getindex(BlockOperator.operators, args...)

#                      ,,,,,,,,,,,,,,,,,,,,,,,,,,,,............................................;[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]])
