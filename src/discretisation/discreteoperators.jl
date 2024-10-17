abstract type DiscreteFractalOperator end

struct DiscreteGalerkinOperator{
    K <: FractalOperator,
    B <: QuasiUniformBasis,
    M <: AbstractMatrix
    } <: DiscreteFractalOperator
    op :: K
    basis :: B
    stiffness_matrix :: M
end

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
            kwargs...)
    Vₕ = construct_p0basis(sio.measure, h_mesh)
    return discretise(sio, Vₕ; kwargs...)
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


function get_galerkin_matrix(sio::AbstractSeparableIntegralOperator{
                        <:AbstractInvariantMeasure, Z},
                    Vₕ::FractalBasis;
                    reps = true
                    ) where Z<:Number
    N = length(Vₕ)

    # initliase Galerkin matrix
    galerkinmatrix = Array{Z}(undef, N, N)

    # get matrix detailing repeated Galerkin entries
    reps ?  galerkinreps = get_galerkinreps(N, sio, Vₕ) :
            galerkinreps = zeros(Int64, N, N)

    # get the information and values of canonical singular integrals
    prepared_singular_vals, singular_indices = getsingularinfo(sio.measure, sio.s, Vₕ.parent_quadrule) 

    # double loop computing inner products, avoiding repeated entries
    @sync for n in 1:N
        @spawn begin
            for m in 1:N
                if @inbounds galerkinreps[m, n] == 0
                    @inbounds galerkinmatrix[m, n] =
                        sesquilinearform(sio, Vₕ[m], Vₕ[n],
                                        prepared_singular_vals, singular_indices)
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

    return galerkinmatrix#DiscreteGalerkinOperator(sio, Vₕ, galerkinmatrix)
end

function get_galerkin_matrix(::IdentityOperator,
                    Vₕ::FractalBasis;
                    # could include 'reps' option here, eventually for homogenous cases
                    )
    return [ϕ⋅ψ for ϕ in Vₕ, ψ in Vₕ]
end

get_galerkin_matrix(S::ScaledOperator, Vₕ::FractalBasis) =
    S.λ * get_galerkin_matrix(S.operator, Vₕ)

get_galerkin_matrix(S::SumOperator, Vₕ::FractalBasis) =
    get_galerkin_matrix(S.operator1, Vₕ) + get_galerkin_matrix(S.operator2, Vₕ)


# generic discretisation function
discretise(K::FractalOperator, Vₕ::FractalBasis; varargs...) =
    DiscreteGalerkinOperator(K, Vₕ, get_galerkin_matrix(K, Vₕ; varargs...))