# define key default parameters

# default constants
DOFS_PER_WAVELENGTH = 20
DOFS_FOR_NONOSCILLATORS = 5

# only makes sense when indices are on homogenous attractor
function vindex_to_scalar(M::Integer, ℓ::Integer, m::AbstractVector{<:Integer})
    n = 1
    for j=1:ℓ
        n += M^(ℓ-j) * (m[j]-1)
    end
    return n
end

abstract type DiscreteFractalOperator end

include("basis.jl")
include("sesquilinearforms.jl")
include("matrixreps.jl")
include("galerkin.jl")
include("collocation.jl")

getdefault_meshwidth(sio::OscillatorySeparableIntegralOperator) =
    2π / abs(DOFS_PER_WAVELENGTH*sio.wavenumber)

getdefault_meshwidth(sio::IntegralOperator) =
    sio.measure.supp.diam / DOFS_FOR_NONOSCILLATORS

# new style which allows us to try different quadrature rules
function discretise(sio::AbstractSingularIntegralOperator;
            h_mesh::Real = getdefault_meshwidth(sio),#sio.measure.supp.diam/5,
            h_quad::Real = 0.0, # quick option for Barycentre rule
            N_quad::Integer = 0,
            quadrule::QuadStruct =
            getdefault_quad_premap(sio.measure, h_mesh, h_quad, N_quad),
            method::Symbol=:galerkin,
            kwargs...)
    Vₕ = construct_quasiuniform_p0basis(sio.measure,
                                        h_mesh,
                                        quadrule)
    if method == :galerkin
        discretise_galerkin(sio, Vₕ; kwargs...)
    elseif method == :collocation
        discretise_collocation(sio, Vₕ, h_quad; kwargs...)
    end
end


struct Projection{B<:FractalBasis, V<:AbstractVector}
    basis::B
    coeffs::V
end

project(Vₕ::FractalBasis, f::Function) = Projection(Vₕ, [ϕₕ(f) for ϕₕ ∈ Vₕ])

# makes sense to rearrange a lot of this, because of how things have evolved
# for example, projections should be inside the corresponding method