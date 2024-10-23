abstract type FractalOperator{M <: Measure} end
include("kernels.jl")
include("integraloperators.jl")

include("operatorpresets.jl")
# include("identityoperators.jl")
include("binaryoperators.jl")