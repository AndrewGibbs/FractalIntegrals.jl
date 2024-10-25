struct DiscreteGalerkinOperator{
    K <: FractalOperator,
    B <: FractalBasis,
    M <: AbstractMatrix
    } <: DiscreteFractalOperator
    op :: K
    basis :: B
    stiffness_matrix :: M
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
discretise_galerkin(K::FractalOperator, Vₕ::FractalBasis; varargs...) =
    DiscreteGalerkinOperator(K, Vₕ, get_galerkin_matrix(K, Vₕ; varargs...))

get_galerkin_matrix(op::SmoothIntegralOperator, dom_basis, codom_basis) =
    [sesquilinearform(op, ϕ,  ψ) for ϕ in dom_basis, ψ in codom_basis]

function get_galerkin_matrix(K::BlockOperator, Vₕ::Tuple{Vararg{FractalBasis}}; varargs...)
    num_bases = length(Vₕ)
    # need to exploit symmetry if possible:
    # need to check we're not computing the transpose block
    block_form = [(m==n ? get_galerkin_matrix(K[m,n], Vₕ[m]; varargs...) :
                        get_galerkin_matrix(K[m,n], Vₕ[m], Vₕ[n]; varargs...))
                        for n in 1:num_bases, m in 1:num_bases]
    full_width = sum(length.(Vₕ))
    # need to intelligently initialise matrix based on types of block matrices
    full_matrix = Matrix{ComplexF64}(undef, full_width, full_width)
    col_loc_start_index = 1
    for n in 1:num_bases # block columns
        # loc_num_rows, loc_num_cols = size(block_form[m,n])
        loc_num_cols = length(Vₕ[n])
        col_inds = col_loc_start_index:(col_loc_start_index+loc_num_cols-1)
        row_loc_start_index = 1
        for m in 1:num_bases # block rows
            loc_num_rows = length(Vₕ[m])
            row_inds = row_loc_start_index:(row_loc_start_index + loc_num_rows - 1)
            full_matrix[row_inds, col_inds] .= block_form[n, m]
            row_loc_start_index = row_inds[end] + 1
        end
        col_loc_start_index = col_inds[end] + 1
    end
    return full_matrix
end

function Base.:\(op::DiscreteGalerkinOperator, f::Function)
    fₕ = project(   op.basis, # basis
                    f)
    coeffs = op.stiffness_matrix \ fₕ.coeffs
    return Projection(op.basis, coeffs)
end