
# split this into a barycentric-centric version and a more general version
struct DiscreteFractalOperator{FO<:FractalOperator,
                                IP<:AbstractInnerProduct,
                                B<:FractalBasis,
                                M<:AbstractMatrix
                                } <: FractalOperator
    op :: FO
    ip :: IP
    basis :: B
    galerkinmatrix :: M
end

function discretise(sio::AbstractSingularIntegralOperator;
                    h_mesh::Number = sio.measure.supp.diam/5,
                    h_quad::Number = h_mesh/5,
                    kwargs...)
    Vₕ = construct_p0basis(sio.measure, h_mesh)
    bip = BarycentreHomogInnerProduct(sio, Vₕ, h_mesh, h_quad)
    return discretise(sio, bip, Vₕ; kwargs...)
end

function count_common_entries(m::AbstractVector{<:Integer}, n::AbstractVector{<:Integer})
    count = 0
    for j in eachindex(m)
        @inbounds m[j] == n[j] ? count +=1 : break
    end
    return count
end

function vindex_to_scalar(M::Integer, ℓ::Integer, m::AbstractVector{<:Integer})
    n = 1
    for j=1:ℓ
        n += M^(ℓ-j) * (m[j]-1)
    end
    return n
end

# function find_redundancies( sio::AbstractSingularIntegralOperator,
#                             Vₕ::FractalBasis)

#     N = length(Vₕ)
#     M = length(sio.measure.supp.ifs)
#     ℓ = round(Int,log(N)/log(M))

#     # this is becoming very general... 
#     # there must be more information we can use about the block structure
#     can_count = 1
#     canonical_index_pairs = Vector{Tuple{Int64, Int64}}(undef,N^2)
#     δ_diff = Vector{eltype(sio.IFS[1].δ)}(undef,N^2)
#     δ_diffs[1] = zero(eltype(sio.IFS[1].δ))
#     δcomps = ..
#     canonical_index_pairs[1] = (1,1)

#     rep_count = 0
#     repeated_index_pairs = Vector{Tuple{Int64 ,Int64}}(undef,N^2)
#     associated_index_pairs = Vector{Tuple{Int64 ,Int64}}(undef,N^2)

#     for (n,ϕ) ∈ enumerate(Vₕ)
#         for (m,ψ) ∈ enumerate(Vₕ)
#             similar = false
#             for j in 1:can_count
#                 if δcomps[m] - δcomps[n] ≈ δ_diffs[j]
#                     similar = true
#                 end
#             end
#             if (L>0) && (m>1 || n>1)
#                 rep_count += 1
#                 # could be done more efficiently, if I can input
#                 # L to vindex_to_scalar and start loop there.
#                 # this would remove allocations and redundant loopands
#                 vn_ = [ones(Int,L); ϕ.vindex[(L+1):end]]
#                 vm_ = [ones(Int,L); ψ.vindex[(L+1):end]]
#                 n_ = vindex_to_scalar(M, ℓ, vn_)
#                 m_ = vindex_to_scalar(M, ℓ, vm_)
#                 repeated_index_pairs[rep_count] = (m,n)
#                 associated_index_pairs[rep_count] = (m_,n_)
#             else
#                 can_count += 1
#                 canonical_index_pairs[can_count] = (m,n)
#             end
#         end
#     end
#     return canonical_index_pairs[1:can_count],
#             repeated_index_pairs[1:rep_count],
#             associated_index_pairs[1:rep_count]
# end

# the code below prevents loops in the repition indices
function repassign!(galerkinreps, fromindex, toindex)
    if galerkinreps[toindex] != 0 # pointing to canonical entry
    #     galerkinreps[fromindex] = toindex
    # else
        while galerkinreps[toindex] != 0 # pointing to another repeated entry
            toindex = galerkinreps[toindex]
        end
    end
    galerkinreps[fromindex] = toindex
end

# this code passes through each 'repeated' entry and checks they point directly to a 'canonical' entry
function cleanreps!(galerkinreps)
    for n in eachindex(galerkinreps)
        if galerkinreps[n] != 0
            repassign!(galerkinreps, n, galerkinreps[n])
        end
    end
end

function selfadjointreps!(galerkinreps, N)
    for n in 1:N
        row_ind = (n-1)*N
        for m in (n+1):N
            repassign!(galerkinreps, row_ind + m, (m-1)*N + n)
            # galerkinreps[m, n] = (m-1)*N + n
        end
    end
end

function fractaldiagblocks!(galerkinreps, N, M, Nglob = N, m_offset = 0, n_offset = 0)
    ℓ = round(Int64, log(N)/log(M))
    for j in 1:ℓ
        block_width = M^(j-1)
        canon_range = 1:block_width
        # first do diagonal blocks
        for M_ in 2:M
            rep_range = (block_width*(M_-1)+1):(block_width*M_)
            for n in canon_range
                for m in canon_range
                    #repassign!(galerkinreps, (rep_range[n]-1)*N + rep_range[m], (n-1)*N + m)
                    galerkinreps[(rep_range[n]-1)*N + rep_range[m]] = (n + n_offset - 1)*Nglob + m + m_offset
                end
            end
        end
    end

    if ℓ>1 # equiv: stop recursing when ℓ==1
        block_width = M^(ℓ-1)
        for m in 1:M
            col_range = (block_width*(m-1)+1):(block_width*m)
            for n in 1:M
                row_range = (block_width*(n-1)+1):(block_width*n)
                if m != n
                    fractaldiagblocks!(@view(galerkinreps[col_range, row_range]),
                                        Int64(N/M), M,
                                        Nglob,
                                        m_offset + col_range[1] - 1,
                                        n_offset + row_range[1] - 1)
                end
            end
        end
    end
end

function fractaloffdiagblocks!(galerkinreps, N, M)
    ℓ = round(Int64, log(N)/log(M))
    # now recursively do off-diagonal blocks
    block_width = M^(ℓ-1)
    for m in 1:M
        col_range = (block_width*(m-1)+1):(block_width*m)
        for n in 1:M
            row_range = (block_width*(n-1)+1):(block_width*n)
            if m != n
                fractaldiagblocks!(@view(galerkinreps[col_range, row_range]), N/M^2, M)
            end
        end
    end
end

function discretise(sio::AbstractSingularIntegralOperator,
                    ip::AbstractInnerProduct,
                    Vₕ::FractalBasis;
                    reps = true
                    )
    N = length(Vₕ)
    M = length(sio.measure.supp.ifs)
    galerkinmatrix = Array{typeof(sio.singularconst)}(undef, N, N)

    # predetermine which entries are repeated in fractal matrix structures
    galerkinreps = zeros(Int64, N, N)
    selfadjointreps!(galerkinreps, N) # not 'fractal' feature - follows from self-adjointness
    if reps
        fractaldiagblocks!(galerkinreps, N, M)
    end
    cleanreps!(galerkinreps)
    
    2+2;
    # first loop avoids repeated entries
    @sync for n in 1:N #@sync 
        @spawn begin #@spawn 
            for m in 1:N
                if @inbounds galerkinreps[m, n] == 0
                    @inbounds galerkinmatrix[m, n] = sesquilinearform(ip, Vₕ[m], Vₕ[n])
                end
            end
        end
    end

    # second loop fills in repeated entries
    @sync for n in 1:N
        @spawn begin
            for m in 1:N
                if @inbounds galerkinreps[m, n] != 0
                    @inbounds galerkinmatrix[m, n] = galerkinmatrix[galerkinreps[m, n]]
                end
            end
        end
    end

    # store matrix as symmetric:
    # most efficient way to invert matrix, also allows to do everything column major

    return DiscreteFractalOperator(sio, ip, Vₕ, Symmetric(galerkinmatrix, :L))
end

struct projection{B<:FractalBasis, V<:AbstractVector}
    basis::B
    coeffs::V
end

function project(bip::BarycentreHomogInnerProduct,
                Vₕ::FractalBasis, 
                f::Function)
    # non-iterate version:
    #Vₕ[n])/Vₕ[n].normalisation for n in eachindex[Vₕ]
    coeffs = [innerproduct(bip, f, ϕₕ)/ϕₕ.normalisation for ϕₕ ∈ Vₕ]
    return projection(Vₕ, coeffs)
end

function Base.:\(op::DiscreteFractalOperator, f::Function)
    fₕ = project(   op.ip, # inner product
                    op.basis, # basis
                    f)
    coeffs = op.galerkinmatrix \ fₕ.coeffs
    return projection(op.basis, coeffs)
end