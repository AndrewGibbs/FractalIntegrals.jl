
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

# following the IFSintegrals style:
# function discretise(sio::AbstractSingularIntegralOperator;
#                     h_mesh::Number = sio.measure.supp.diam/5,
#                     h_quad::Number = h_mesh/5,
#                     kwargs...)
#     Vₕ = construct_p0basis(sio.measure, h_mesh)
#     quadpts, quadweights = barycentre_quadrule(sio.measure, h_quad)
#     ip = InnerProduct(sio, Vₕ, quadpts, quadweights)
#     return discretise(sio, ip, Vₕ; kwargs...)
# end

getdefault_meshwidth(sio::OscillatorySingularIntegralOperator) = 2π / abs(10*sio.wavenumber)
getdefault_meshwidth(sio::SingularIntegralOperator) = sio.diam / 5
function getdefault_quad(sio::AbstractSingularIntegralOperator; h_quad::Real = 0.0)
    if h_quad > 0
        barycentre_quadrule(sio.measure, h_quad)
    elseif sio.measure.supp.n == 1
        # replace this with Gauss rule when I've coded it
        @warn("Need to replace this option with Mantica Gauss")
        barycentre_quadrule(sio.measure, getdefault_meshwidth(sio))
    else
        barycentre_quadrule(sio.measure, getdefault_meshwidth(sio))
    end
end

# new style which allows us to try different quadrature rules
function discretise(sio::AbstractSingularIntegralOperator;
                    h_mesh::Real = getdefault_meshwidth(sio),#sio.measure.supp.diam/5,
                    h_quad::Real = 0.0, # quick option for Barycentre rule
                    quadrule::Tuple{AbstractVector,AbstractVector} = getdefault_quad(sio, h_quad = h_quad),
                    kwargs...)
    Vₕ = construct_p0basis(sio.measure, h_mesh)
    # quadpts, quadweights = barycentre_quadrule(sio.measure, h_quad)
    ip = InnerProduct(sio, Vₕ, quadrule[1], quadrule[2])
    return discretise(sio, ip, Vₕ; kwargs...)
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

function get_galerkinreps(N::Integer, sio::AbstractSingularIntegralOperator)
        M = length(sio.measure.supp.ifs)
        # predetermine which entries are repeated in fractal matrix structures
        galerkinreps = zeros(Int64, N, N)
        if sio.selfadjoint
            selfadjointreps!(galerkinreps, N) # not 'fractal' feature - follows from self-adjointness
        end
        if isa(sio.measure.supp.ifs, AbstractVector{<:TranslatingSimilarity})
            fractaldiagblocks!(galerkinreps, N, M)
        end
        cleanreps!(galerkinreps)
        return galerkinreps
end

function discretise(sio::AbstractSingularIntegralOperator{
                            <:AbstractInvariantMeasure{
                                <:AbstractAttractor{T, R}}, 
                            Z},
                    ip::AbstractInnerProduct,
                    Vₕ::FractalBasis;
                    reps = true
                    ) where {
                    T, R<:Real, Z<:Number}
    N = length(Vₕ)
    μ = sio.measure
    
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

struct Projection{B<:FractalBasis, V<:AbstractVector}
    basis::B
    coeffs::V
end

function project(bip::InnerProduct,
                Vₕ::FractalBasis, 
                f::Function)
    # non-iterate version:
    #Vₕ[n])/Vₕ[n].normalisation for n in eachindex[Vₕ]
    coeffs = [innerproduct(bip, f, ϕₕ)/ϕₕ.normalisation for ϕₕ ∈ Vₕ]
    return Projection(Vₕ, coeffs)
end

function Base.:\(op::DiscreteFractalOperator, f::Function)
    fₕ = project(   op.ip, # inner product
                    op.basis, # basis
                    f)
    coeffs = op.galerkinmatrix \ fₕ.coeffs
    return Projection(op.basis, coeffs)
end