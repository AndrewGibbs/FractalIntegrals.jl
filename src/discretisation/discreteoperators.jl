# split this into a barycentric-centric version and a more general version
struct DiscreteFractalOperator{
        FO<:FractalOperator,
        IP<:AbstractInnerProduct,
        B<:FractalBasis,
        M<:AbstractMatrix
        } <: FractalOperator
    op :: FO
    ip :: IP
    basis :: B
    galerkinmatrix :: M
end

getdefault_meshwidth(sio::OscillatorySingularIntegralOperator) =
    2π / abs(DOFS_PER_WAVELENGTH*sio.wavenumber)

getdefault_meshwidth(sio::IntegralOperator) =
    sio.measure.supp.diam / DOFS_FOR_NONOSCILLATORS


# getdefault_quadwidth(Γ::AbstractAttractor) =
#     maximum(sₘ.ρ for sₘ in Γ.ifs)^QUAD_EXTRA_LEVELS

# getdefault_quadwidth(sio::IntegralOperator) =
#     # getdefault_meshwidth(sio) *
#         diam(sio.measure) * maximum(sₘ.ρ for sₘ in sio.measure.supp.ifs)^QUAD_EXTRA_LEVELS

# getdefault_gaussorder(::IntegralOperator) = QUAD_DEFAULT_GAUSS

# function getdefault_quad(sio::AbstractSingularIntegralOperator,
#                         h_mesh = diam(sio.measure);
#                         h_quad::Real = 0.0,
#                         N_quad::Integer = 0)
#     if h_quad > 0
#         x, w = barycentre_quadrule(sio.measure, h_quad*diam(sio.measure)/h_mesh)
#     elseif sio.measure.supp.n == 1
#         if N_quad >= 1
#             x, w = gauss_quadrule(sio.measure, N_quad)
#         else
#             x, w = gauss_quadrule(sio.measure, getdefault_gaussorder(sio))
#         end
#     else
#         x, w = barycentre_quadrule(sio.measure, getdefault_quadwidth(sio))
#     end
#     return x, w
# end

getdefault_quad_premap(μ, h_mesh, h_quad, N_quad) =
    getdefault_quad(μ, h_quad*diam(μ)/h_mesh, N_quad)


function getdefault_quad(μ::AbstractInvariantMeasure,
                        h_quad,
                        N_quad)
        if h_quad > 0
            x, w = barycentre_quadrule(μ, h_quad)
        elseif μ.supp.n == 1
            if N_quad >= 1
                x, w = gauss_quadrule(μ, N_quad)
            else
                x, w = gauss_quadrule(μ, QUAD_DEFAULT_GAUSS)
            end
        else
            x, w = barycentre_quadrule(μ, default_barywidth(μ))
        end
    return QuadStruct(x, w)
end

# function getdefault_quad(sio::AbstractSingularIntegralOperator;
#                         h_quad::Real = 0.0,
#                         N_quad::Integer = 0)
#     if h_quad > 0
#         x, w = barycentre_quadrule(sio.measure, h_quad)
#     elseif sio.measure.supp.n == 1
#         if N_quad >= 1
#             x, w = gauss_quadrule(sio.measure, N_quad)
#         else
#             x, w = gauss_quadrule(sio.measure, getdefault_gaussorder(sio))
#         end
#     else
#         x, w = barycentre_quadrule(sio.measure, getdefault_quadwidth(sio))
#     end
#     return x, w
# end

# new style which allows us to try different quadrature rules
function discretise(sio::AbstractSingularIntegralOperator;
    h_mesh::Real = getdefault_meshwidth(sio),#sio.measure.supp.diam/5,
    h_quad::Real = 0.0, # quick option for Barycentre rule
    N_quad::Integer = 0,
    quadrule::QuadStruct =
            getdefault_quad_premap(sio.measure, h_mesh, h_quad, N_quad),
            kwargs...)
    Vₕ = construct_p0basis(sio.measure, h_mesh)
    # quadpts, quadweights = barycentre_quadrule(sio.measure, h_quad)
    ip = InnerProduct(sio, Vₕ, quadrule)
    return discretise(sio, ip, Vₕ; kwargs...)
end

function convert_quad_to_tuple(Q::QuadStruct)
    (; nodes, weights) = Q
    return (nodes, weights)
end

function count_common_entries(m::AbstractVector{<:Integer}, n::AbstractVector{<:Integer})
    count = 0
    for j in eachindex(m)
        @inbounds m[j] == n[j] ? count +=1 : break
    end
    return count
end

function discretise(sio::AbstractSingularIntegralOperator{
                        <:AbstractInvariantMeasure{
                            <:AbstractAttractor{R, T}}, 
                        Z},
                    ip::AbstractInnerProduct,
                    Vₕ::FractalBasis;
                    reps = true
                    ) where {
                    T, R<:Real, Z<:Number}
N = length(Vₕ)

# get element type of Galerkin matrix
ElementType = promote_type(Z, eltype(T), R)

# initliase Galerkin matrix
galerkinmatrix = Array{ElementType}(undef, N, N)

# get matrix detailing repeated Galerkin entries
reps ? galerkinreps =
    get_galerkinreps(N, sio) :
    galerkinreps = zeros(Int64, N, N)

# double loop computing inner products, avoiding repeated entries
@sync for n in 1:N #@sync 
    @spawn begin #@spawn 
        for m in 1:N
            if @inbounds galerkinreps[m, n] == 0
                @inbounds galerkinmatrix[m, n] = sesquilinearform(ip, Vₕ[m], Vₕ[n])
            end
        end
    end
end

# second loop filling in repeated entries
@sync for n in 1:N
    @spawn begin
        for m in 1:N
            if @inbounds galerkinreps[m, n] != 0
                @inbounds galerkinmatrix[m, n] = galerkinmatrix[galerkinreps[m, n]]
            end
        end
    end
end

    return DiscreteFractalOperator(sio, ip, Vₕ, galerkinmatrix)
end