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
                        <:AbstractInvariantMeasure, Z},
                    ip::AbstractInnerProduct,
                    Vₕ::FractalBasis;
                    reps = true
                    ) where Z
    N = length(Vₕ)

    # initliase Galerkin matrix
    galerkinmatrix = Array{Z}(undef, N, N)

    # get matrix detailing repeated Galerkin entries
    reps ?  galerkinreps = get_galerkinreps(N, sio, Vₕ) :
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