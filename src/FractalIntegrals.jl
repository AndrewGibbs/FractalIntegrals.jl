module FractalIntegrals

import StaticArrays: SVector, SMatrix, mul!
import LinearAlgebra: norm, det, dot
import LinearAlgebra.I as IdMat
import SpecialFunctions: hankelh1
export Similarity

include("affinemaps.jl")
include("fractalmeasures.jl")
include("barycentrerule.jl")
include("fastkernels.jl")
include("senergy.jl")
include("basis.jl")
include("operators.jl")
include("operatorpresets.jl")
include("dualpairings.jl")

end
