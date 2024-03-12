module FractalIntegrals

import StaticArrays: SVector, SMatrix, mul!
import LinearAlgebra: norm, det
import LinearAlgebra.I as IdMat
export Similarity

include("affinemaps.jl")
include("fractalmeasures.jl")
include("barycentrerule.jl")
# include("senergy.jl")

end
