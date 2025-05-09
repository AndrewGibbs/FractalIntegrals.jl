# define key default parameters

# default constants
DOFS_PER_WAVELENGTH = 20
DOFS_FOR_NONOSCILLATORS = 5

# only makes sense when indices are on homogenous attractor
function vindex_to_scalar(M::Integer, ℓ::Integer, m::AbstractVector{<:Integer})
    n = 1
    for j = 1:ℓ
        n += M^(ℓ - j) * (m[j] - 1)
    end
    return n
end

abstract type DiscreteFractalOperator end

include("basis.jl")
include("sesquilinearforms.jl")
include("matrixreps.jl")
include("galerkin.jl")
include("collocation.jl")

function getdefault_meshwidth(sio::OscillatorySeparableIntegralOperator)
    return 2π / abs(DOFS_PER_WAVELENGTH * sio.wavenumber)
end

function getdefault_meshwidth(sio::IntegralOperator)
    return sio.measure.supp.diam / DOFS_FOR_NONOSCILLATORS
end

# new style which allows us to try different quadrature rules
function discretise(
    sio::FractalOperator;
    h_mesh::Real = getdefault_meshwidth(sio),#sio.measure.supp.diam/5,
    h_quad::Real = 0.0, # quick option for Barycentre rule
    N_quad::Integer = 0,
    quadrule::QuadStruct = getdefault_quad_premap(sio.measure, h_mesh, h_quad, N_quad),
    method::Symbol = :galerkin,
    kwargs...,
)
    Vₕ = construct_quasiuniform_p0basis(sio.measure, h_mesh, quadrule)
    if method == :galerkin
        discretise_galerkin(sio, Vₕ; kwargs...)
    elseif method == :collocation
        discretise_collocation(sio, Vₕ, h_quad; kwargs...)
    else
        error("method not recognised")
    end
end

struct Projection{B<:FractalBasis,V<:AbstractVector}
    basis::B
    coeffs::V
end

project(Vₕ::FractalBasis, f::Function) = Projection(Vₕ, [ϕₕ(f) for ϕₕ in Vₕ])

# (p::Projection)(x) = sum(p.basis[n](x)*p.coeffs[n] for n=1:length(p.basis))

# makes sense to rearrange a lot of this, because of how things have evolved
# for example, projections should be inside the corresponding method
