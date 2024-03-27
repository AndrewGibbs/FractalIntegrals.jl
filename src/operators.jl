abstract type FractalOperator end
abstract type AbstractSingularIntegralOperator <: FractalOperator end


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