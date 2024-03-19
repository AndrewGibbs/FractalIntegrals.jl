abstract type FractalOperator end
abstract type SingularIntegralOperator <: FractalOperator end

struct SeparableIntegralOperator{M <: AbstractInvariantMeasure,
                                S <: Number,
                                Z <: Number
                                } <: SingularIntegralOperator
    measure::M
    kernel::Function
    lipschitzpart::Function
    s::S
    singularconst::Z
end