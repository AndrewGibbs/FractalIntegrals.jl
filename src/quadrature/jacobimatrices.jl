function map_coeffs_to_tilde_coeffs(coeffs::Vector{T}, rₙ::T) where {T<:Real}
    n = length(coeffs)
    tilde_coeffs = zeros(T, n)
    tilde_coeffs[1:(n-1)] = coeffs[1:(n-1)] .* rₙ
    tilde_coeffs[n] = coeffs[n]
    return tilde_coeffs
end

function map_tilde_coeffs_to_coeffs(tilde_coeffs::Vector{T}, rₙ::T) where {T<:Real}
    n = length(tilde_coeffs)
    coeffs = zeros(T, n)
    coeffs[1:(n-1)] = tilde_coeffs[1:(n-1)] ./ rₙ
    coeffs[n] = tilde_coeffs[n]
    return coeffs
end

# now overload both of the above functions with full i=1,…,M inputs
function map_coeffs_to_tilde_coeffs(coeffs::Vector{Vector{T}}, rₙ::T) where {T<:Real}
    M = length(coeffs)
    tilde_coeffs = [map_coeffs_to_tilde_coeffs(coeffs[i], rₙ) for i = 1:M]
    return tilde_coeffs
end
function map_tilde_coeffs_to_coeffs(tilde_coeffs::Vector{Vector{T}}, rₙ::T) where {T<:Real}
    M = length(tilde_coeffs)
    coeffs = [map_tilde_coeffs_to_coeffs(tilde_coeffs[i], rₙ) for i = 1:M]
    return coeffs
end

function get_next_modified_coeffs_from_coeffs(
    Γⁿ⁻¹::AbstractVector{T},
    Γⁿ⁻²::AbstractVector{T},
    A::AbstractVector{T},
    r::AbstractVector{T},
    sₘ::AbstractSimilarity,
) where {T<:Real}
    # rule of thumb: square bracket indices should always be incremented by one compared with what's written above,
    # to account for the zero index

    # Conflicting notation between Mantica '96 and myself, so this step should switch everything to his notation:
    δᵢ = sₘ.ρ
    βᵢ = sₘ.δ
    # that's better.

    n = length(Γⁿ⁻¹) # one more because of the zeroth entry, one less because we're incrementing
    modΓⁿ = zeros(T, n + 1)

    #  modΓⁿ⁻¹ = map_coeffs_to_tilde_coeffs(Γⁿ⁻¹,r[end])
    #  modΓⁿ⁻² = map_coeffs_to_tilde_coeffs(Γⁿ⁻²,r[end])

    function γ(ℓ::Integer, j::Integer)
        if ℓ == j
            return A[ℓ+1] * Γⁿ⁻¹[ℓ+1]
        elseif ℓ == j - 1
            return r[j+1] * Γⁿ⁻¹[ℓ+1]
        elseif j == ℓ - 1
            return r[ℓ+1] * Γⁿ⁻¹[ℓ+1]
        else
            return 0.0
        end
    end

    function get_trickier_sum(j)
        trickier_sum = 0.0
        for ℓ = max(0, j - 1):min(n - 1, j + 1) # limit this to the non-zero stuff
            trickier_sum += γ(ℓ, j)
        end
        return trickier_sum
    end

    for j = 0:(n-2)
        modΓⁿ[j+1] = (βᵢ - A[n]) * Γⁿ⁻¹[j+1] + δᵢ * get_trickier_sum(j) - r[n] * Γⁿ⁻²[j+1]
    end

    # now deal with j = n-1
    modΓⁿ[n] = (βᵢ - A[n]) * Γⁿ⁻¹[n] + δᵢ * get_trickier_sum(n - 1)

    # now deal with j = n, which is not in orthanormal function
    modΓⁿ[n+1] = δᵢ * Γⁿ⁻¹[n]

    return modΓⁿ
end

function get_rₙ(
    modΓⁿ::Vector{Vector{T}},
    Γⁿ⁻¹::Vector{Vector{T}},
    A::Vector{T},
    r::Vector{T},
    μ::AbstractInvariantMeasure,
) where {T<:Real}
    n = length(modΓⁿ[1]) - 1
    M = length(μ.supp.ifs)
    B = zeros(T, M)
    C = zeros(T, M)
    #  D = zeros(Float64,M)
    D_nominator = zeros(T, M)

    for i = 1:M
        for ℓ = 0:(n-1)
            B[i] +=
                (μ.supp.ifs[i].δ + μ.supp.ifs[i].ρ * A[ℓ+1]) * modΓⁿ[i][ℓ+1] * Γⁿ⁻¹[i][ℓ+1]
        end
        for ℓ = 0:(n-2) # Not used for n=1
            C[i] += r[ℓ+2] * (modΓⁿ[i][ℓ+1] * Γⁿ⁻¹[i][ℓ+2] + modΓⁿ[i][ℓ+2] * Γⁿ⁻¹[i][ℓ+1])
        end
        C[i] *= μ.supp.ifs[i].ρ
        D_nominator[i] = μ.supp.ifs[i].ρ * modΓⁿ[i][n+1] * Γⁿ⁻¹[i][n]
    end

    rₙ² = (μ.weights' * (B + C)) / (1 - (μ.weights' * D_nominator))

    return sqrt(max(rₙ², zero(T)))
end

function get_Aₙ(
    coeffs_this_level::Vector{Vector{T}},
    A::Vector{T},
    r::Vector{T},
    μ::AbstractInvariantMeasure,
) where {T<:Real}
    M = length(μ.supp.ifs)
    n = length(coeffs_this_level[1]) - 1
    big_sum = zeros(M)
    denominator = 0.0
    for i = 1:M
        big_sum[i] += coeffs_this_level[i][n+1]^2 * μ.supp.ifs[i].δ
        for m = 0:(n-1)
            big_sum[i] +=
                coeffs_this_level[i][m+1]^2 * (μ.supp.ifs[i].δ + μ.supp.ifs[i].ρ * A[m+1])# +...
            big_sum[i] +=
                μ.supp.ifs[i].ρ *
                (2 * coeffs_this_level[i][m+1] * coeffs_this_level[i][m+2]) *
                r[m+2]
        end
        denominator += μ.weights[i] * coeffs_this_level[i][n+1]^2 * μ.supp.ifs[i].ρ
    end
    return (μ.weights' * big_sum) / (1 - denominator)
end

function getjacobimatrix(
    μ::AbstractInvariantMeasure{1,M,T,<:AbstractAttractor},
    N::Integer,
) where {M,T<:Real}

    # initialisation
    A = zeros(T, N + 1)
    r = zeros(T, N + 1)
    J = zeros(T, N + 1, N + 1)
    # M = length(μ.supp.ifs)
    coeffs_one_below = [[one(T)] for _ = 1:M]
    coeffs_two_below = [[zero(T)] for _ = 1:M]
    A[1] = μ.suppmeasure * get_barycentre(μ) # checks out, given def'n of barycentre, should be ∫_Γ x dμ(x)

    # iteration
    for n = 1:N #we've done n=0 above
        # initialise this level
        coeffs_this_level = [zeros(T, n + 1) for _ = 1:M]

        # step one
        modified_coeffs = [zeros(T, n + 1) for _ = 1:M]
        for m = 1:M
            modified_coeffs[m] = get_next_modified_coeffs_from_coeffs(
                coeffs_one_below[m],
                coeffs_two_below[m],
                A,
                r,
                μ.supp.ifs[m],
            )
        end

        # step two
        r[n+1] = get_rₙ(modified_coeffs, coeffs_one_below, A, r, μ)

        # step three
        coeffs_this_level = map_tilde_coeffs_to_coeffs(modified_coeffs, r[n+1])

        # step four
        A[n+1] = get_Aₙ(coeffs_this_level, A, r, μ)

        coeffs_two_below = coeffs_one_below
        coeffs_one_below = coeffs_this_level
    end

    # now make the Jacobi matrix
    for i = 0:(N-1)
        J[i+1, i+1] = A[i+1]
        J[i+1, i+2] = r[i+2]
        J[i+2, i+1] = r[i+2]
    end
    J[N+1, N+1] = A[N+1]

    return J
end
