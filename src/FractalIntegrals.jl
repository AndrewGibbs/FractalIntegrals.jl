module FractalIntegrals

import StaticArrays: SVector, SMatrix, mul!
import LinearAlgebra: norm, det, dot, Symmetric, kron, eigvecs, eigvals
import LinearAlgebra.I as IdMat
import SpecialFunctions: hankelh1
import Roots: find_zero, Bisection
import Base.Threads: @spawn, nthreads
export Similarity, FractaPresets
import ChunkSplitters: chunks
using Plots

export getfractal

include("affinemaps.jl")
include("fractalmeasures.jl")
include("boundingball.jl")
include("barycentrerule.jl")
include("fastkernels.jl")
include("senergy.jl")
include("basis.jl")
include("operators.jl")
include("operatorpresets.jl")
include("innerproducts.jl")
include("projections.jl")
# export discretise
include("fractalpresets.jl")
# export cantorset, cantordust
include("potentials.jl")
include("plots.jl")
include("jacobimatrices.jl")

end
