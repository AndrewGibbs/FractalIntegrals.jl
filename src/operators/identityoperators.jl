abstract type FractalOperator end
struct IdentityOperator end
struct ScaledIdendityOperator{T} <: FractalOperator
    λ::T
end
