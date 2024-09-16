abstract type FractalOperator end
abstract type IntegralOperator{M<:AbstractInvariantMeasure, Z} <: FractalOperator end
abstract type AbstractSingularIntegralOperator{M, Z} <: IntegralOperator{M, Z} end

include("kernels.jl")

struct SingularIntegralOperator{
        M <: AbstractInvariantMeasure,
        F1 <: Function,
        F2 <: Function,
        S <: Real,
        Z <: Number
        } <: AbstractSingularIntegralOperator{M, Z}
    measure::M
    kernel::F1
    lipschitzpart::F2
    s::S
    singularconst::Z
    symmetric::Bool   
end

struct OscillatorySingularIntegralOperator{
        M <: AbstractInvariantMeasure,
        F1 <: Function,
        F2 <: Function,
        S <: Real,
        Z <: Number,
        } <: AbstractSingularIntegralOperator{M, Z}
    measure::M
    kernel::F1
    lipschitzpart::F2
    s::S
    singularconst::Z
    symmetric::Bool
    wavenumber::S
end

SingularIntegralOperator(Γ::AbstractAttractor, vargs...) = SingularIntegralOperator(HausdorffMeasure(Γ), vargs...)
OscillatorySingularIntegralOperator(Γ::AbstractAttractor, vargs...) = OscillatorySingularIntegralOperator(Γ::AbstractAttractor, vargs...)

include("operatorpresets.jl")