abstract type FractalOperator end
abstract type IntegralOperator <: FractalOperator end
abstract type AbstractSingularIntegralOperator <: IntegralOperator end


struct SingularIntegralOperator{M <: AbstractInvariantMeasure,
                                F1 <: Function,
                                F2 <: Function,
                                S <: Number,
                                Z <: Number
                                } <: AbstractSingularIntegralOperator
    measure::M
    kernel::F1
    lipschitzpart::F2
    s::S
    singularconst::Z
end

struct SelfAdjointSingularIntegralOperator{M <: AbstractInvariantMeasure,
                                            F1 <: Function,
                                            F2 <: Function,
                                            S <: Number,
                                            Z <: Number
                                            } <: AbstractSingularIntegralOperator
    measure::M
    kernel::F1
    lipschitzpart::F2
    s::S
    singularconst::Z
end

SingularIntegralOperator(Γ::AbstractAttractor, vargs...) = SingularIntegralOperator(HausdorffMeasure(Γ), vargs...)