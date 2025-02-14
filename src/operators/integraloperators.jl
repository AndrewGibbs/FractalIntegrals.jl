abstract type IntegralOperator{M<:Measure,Z<:Number} <: FractalOperator{M} end
abstract type AbstractSingularIntegralOperator{M<:AbstractInvariantMeasure,Z} <:
              IntegralOperator{M,Z} end
abstract type AbstractSeparableIntegralOperator{M,Z} <:
              AbstractSingularIntegralOperator{M,Z} end

struct SeparableIntegralOperator{
    M<:AbstractInvariantMeasure,
    F1<:Function,
    F2<:Function,
    S<:Real,
    Z<:Number,
} <: AbstractSeparableIntegralOperator{M,Z}
    measure::M
    kernel::F1
    lipschitzpart::F2
    s::S
    singularconst::Z
    symmetric::Bool
end

struct SmoothIntegralOperator{
    Mdom<:AbstractInvariantMeasure,
    Mcodom<:AbstractInvariantMeasure,
    F<:Function,
    Z<:Number,
} <: IntegralOperator{Mdom,Z}
    measure::Mdom
    measure_codom::Mcodom
    kernel::F
    symmetric::Bool
end

struct SingularIntegralOperator{M<:AbstractInvariantMeasure,F<:Function,Z<:Number} <:
       AbstractSingularIntegralOperator{M,Z}
    measure::M
    kernel::F
    symmetric::Bool
end

# ultimately, should bin off the below type, its too messy
struct OscillatorySeparableIntegralOperator{
    M<:AbstractInvariantMeasure,
    F1<:Function,
    F2<:Function,
    S<:Real,
    Z<:Number,
} <: AbstractSeparableIntegralOperator{M,Z}
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
