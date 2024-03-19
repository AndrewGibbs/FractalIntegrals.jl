module FractalIntegrals

import StaticArrays: SVector, SMatrix, mul!
import LinearAlgebra: norm, det
import LinearAlgebra.I as IdMat
export Similarity

include("affinemaps.jl")
include("fractalmeasures.jl")
include("barycentrerule.jl")
include("fastkernels.jl")
include("senergy.jl")
include("basis.jl")

end
