
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