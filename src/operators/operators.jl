abstract type FractalOperator{M <: AbstractInvariantMeasure} end
abstract type IntegralOperator{M, Z} <: FractalOperator{M} end
abstract type AbstractSingularIntegralOperator{M, Z} <: IntegralOperator{M, Z} end
abstract type AbstractSeparableIntegralOperator{M, Z} <: AbstractSingularIntegralOperator{M, Z} end

include("kernels.jl")

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

include("operatorpresets.jl")
# include("identityoperators.jl")
include("binaryoperators.jl")