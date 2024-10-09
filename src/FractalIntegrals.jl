module FractalIntegrals

import StaticArrays: SVector, SMatrix, mul!
import LinearAlgebra: norm, det, dot, Symmetric, kron, eigvecs, eigvals
import LinearAlgebra.I as IdMat
import SpecialFunctions: hankelh1, gamma as gammafn
import Roots: find_zero, Bisection
import Base.Threads: @spawn, nthreads
import ChunkSplitters: chunks
using Plots#, Documenter

export Similarity, Attractor, HausdorffMeasure, InvariantMeasure,
        getfractal, FractaPresets

include("fractalgeometry/fractalgeometry.jl")
include("vectorindices/vectorindices.jl")
include("quadrature/quadrature.jl")
include("integraloperators/integraloperators.jl")
include("discretisation/discretisation.jl")

# # less 'core' routines
include("potentials/potentials.jl")
include("fractalpresets/fractalpresets.jl")
include("plots/plots.jl")

end
