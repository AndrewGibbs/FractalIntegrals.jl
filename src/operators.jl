abstract type FractalOperator end
abstract type AbstractSingularIntegralOperator <: FractalOperator end


struct SelfAdjointSingularIntegralOperator{M <: AbstractInvariantMeasure,
                                S <: Number,
                                Z <: Number
                                } <: AbstractSingularIntegralOperator
    measure::M
    kernel::Function
    lipschitzpart::Function
    s::S
    singularconst::Z
end