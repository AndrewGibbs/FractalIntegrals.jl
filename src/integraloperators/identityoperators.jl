abstract type FractalOperator end
struct ScaledIdendityOperator{T} <: FractalOperator
    λ::T
end
