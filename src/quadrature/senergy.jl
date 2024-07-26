# The below variation of vcat is needed for our convention that Γ₀:=Γ
vcat_(x,y) = vcat(x[x .!= 0],y[y .!= 0])
hcat_(x,y) = hcat(x[x .!= 0],y[y .!= 0])

# to be similar, two similarities must pass three tests:
similar_test_one(sₘ::AbstractSimilarity,
                sₙ::AbstractSimilarity,
                sₘ_::AbstractSimilarity,
                sₙ_::AbstractSimilarity
                ) = isapprox(sₘ.ρ/sₙ.ρ, sₘ_.ρ/sₙ_.ρ, atol=100*eps())

similar_test_two(sₘ::AbstractSimilarity,
                sₙ::AbstractSimilarity,
                sₘ_::AbstractSimilarity,
                sₙ_::AbstractSimilarity,
                T::AbstractSimilarity,
                T_::AbstractSimilarity
                ) = isapprox(sₘ.A*T.A*inv(sₙ.A) , sₘ_.A*T_.A*inv(sₙ_.A), atol=100*eps())

similar_test_three(sₘ::AbstractSimilarity,
                sₙ::AbstractSimilarity,
                sₘ_::AbstractSimilarity,
                sₙ_::AbstractSimilarity,
                T::AbstractSimilarity,
                T_::AbstractSimilarity
                )=isapprox(sₘ.δ .- sₘ_.δ .- sₘ.ρ*sₘ.A*(T.A/sₙ.ρ*inv(sₙ.A)*sₙ.δ .- T.δ), .- sₘ_.ρ*sₘ_.A*(T_.A/sₙ_.ρ*inv(sₙ_.A)*sₙ_.δ .- T_.δ), atol=100*eps())

# now account for the special case of no rotation
# similar_test_two(sₘ::TranslatingSimilarity,
#                 sₙ::TranslatingSimilarity,
#                 sₘ_::TranslatingSimilarity,
#                 sₙ_::TranslatingSimilarity,
#                 T::AbstractInvariantMap,
#                 T_::AbstractInvariantMap
#                 ) = true

# similar_test_three(sₘ::TranslatingSimilarity,
#                 sₙ::TranslatingSimilarity,
#                 sₘ_::TranslatingSimilarity,
#                 sₙ_::TranslatingSimilarity,
#                 T::AbstractInvariantMap,
#                 T_::AbstractInvariantMap
#                 ) = isapprox(sₘ.δ .- sₘ_.δ .- sₘ.ρ*T.A*(sₙ.δ .- sₙ_.δ), atol=100*eps())

function check_if_similar(  Γ::AbstractAttractor{R, V},
                            m::AbstractVector{<:Integer},
                            n::AbstractVector{<:Integer},
                            m_::AbstractVector{<:Integer},
                            n_::AbstractVector{<:Integer},
                            G::AbstractVector{<:AbstractSimilarity},
                            G_::AbstractVector{<:AbstractSimilarity}
                            ) where {R, V}
    # get shorthand for IFS
    S = Γ.ifs
    test_one_pass = false
    pass_all = false
    ρ = zero(S[1].ρ) # initialise

    # define identity similarity, which is a workaround for index [0]
    s₀ = one(Γ.ifs[1])

    m != [0] ? sₘ = simmulticomp(S,m) : sₘ = s₀
    n != [0] ? sₙ = simmulticomp(S,n) : sₙ = s₀
    m_ !=[0] ? sₘ_ = simmulticomp(S,m_) : sₘ_ = s₀
    n_ !=[0] ? sₙ_ = simmulticomp(S,n_) : sₙ_ = s₀

    if similar_test_one(sₘ, sₙ, sₘ_, sₙ_)
        ρ = sₘ.ρ/sₙ.ρ
        test_one_pass = true
    end

    if test_one_pass
        for T ∈ G, T_ ∈ G_
            similar_test_two(sₘ, sₙ, sₘ_, sₙ_, T, T_) && similar_test_three(sₘ, sₙ, sₘ_, sₙ_, T, T_) ? pass_all = true : pass_all = false

            # if tests are satisfies for some particular group element pair, can break loop early:
            if pass_all
                break
            end
        end
    end
    return pass_all, ρ
end

# HAVE EDITED DOWN TO HERE

function check_for_similar_integrals(Γ::AbstractAttractor,
                                    X::Vector{<:Tuple{Vector{<:Integer}, Vector{<:Integer}}}, 
                                    mcat::Vector{<:Integer},
                                    mcat_::Vector{<:Integer},
                                    G₁::AbstractVector{<:AbstractSimilarity},
                                    G₂::AbstractVector{<:AbstractSimilarity},
                                    fubini_flag::Bool
                                    )
    is_X_similar = false
    similar_index = -1
    proportionality_const = zero(Γ.ifs[1].ρ)

    # should compactly write the following as a function, it's almost repeated
    for ∫∫_index in 1:length(X)
        ∫∫_indices_higher_level = X[∫∫_index]
        this_is_X_similar = true
        this_is_X_similar, ρ = check_if_similar(Γ, ∫∫_indices_higher_level[1], mcat, ∫∫_indices_higher_level[2], mcat_, G₁, G₂)
        
        if !this_is_X_similar && fubini_flag
            this_is_X_similar, ρ = check_if_similar(Γ, ∫∫_indices_higher_level[2], mcat, ∫∫_indices_higher_level[1], mcat_, G₁, G₂)
        end

        # if we've found a similarity, terminate the process early
        if this_is_X_similar
            is_X_similar = true
            proportionality_const = ρ
            similar_index = ∫∫_index
            break
        end
    end
    return is_X_similar, proportionality_const, similar_index
end

function convert_vector_index_to_integer_index(m::AbstractVector{<:Integer}, M::Integer)
    ℓ = length(m)
    m_integer = 0
    for ℓ_=1:(ℓ-1)
        m_integer += M^(ℓ-ℓ_)*(m[ℓ_]-1)
    end
    m_integer += m[end]
    return m_integer
end

function check_for_ℓ_singular_integrals(Γ::AbstractAttractor, m::Vector{<:Integer}, n::Vector{<:Integer})
    is_singular = -1

    if m==n || m==[0] || n==[0]
        is_singular = 1
    else
        # get important bits
        M = length(Γ.ifs)
        Γ_singularities = Γ.connectedness
        ℓ_depth = Int64(round(log(size(Γ_singularities)[1])/log(M)))

        m_ℓ_remainder_depth = ℓ_depth-length(m)
        n_ℓ_remainder_depth = ℓ_depth-length(n)
        if m_ℓ_remainder_depth>=0 && n_ℓ_remainder_depth>=0
            m_start = convert_vector_index_to_integer_index([m; ones(Int64,m_ℓ_remainder_depth)], M::Int64)
            m_end = convert_vector_index_to_integer_index([m; M*ones(Int64,m_ℓ_remainder_depth)], M::Int64)
            m_range = m_start:m_end
            
            n_start = convert_vector_index_to_integer_index([n; ones(Int64,n_ℓ_remainder_depth)], M::Int64)
            n_end = convert_vector_index_to_integer_index([n; M*ones(Int64,n_ℓ_remainder_depth)], M::Int64)
            n_range = n_start:n_end

            sum(Γ_singularities[m_range,n_range])>0 ? is_singular = 1 : is_singular = 0
        end
    end

    return is_singular
end

function compose_weights(pw::AbstractVector{<:Real}, m::AbstractVector{<:Integer})
    m  != [0] ? pₘ = prod(pw[m]) : pₘ = one(eltype(pw))
    return pₘ
end

function similar_scaler(ρ::Real,
                        s::Real,
                        m::AbstractVector{<:Integer},
                        m_::AbstractVector{<:Integer},
                        n::AbstractVector{<:Integer},
                        n_::AbstractVector{<:Integer},
                        pw₁::AbstractVector{<:Real},
                        pw₂::AbstractVector{<:Real})

    # account for convention Γ₀:=Γ
    pₘ = compose_weights(pw₁, m)
    pₘ_ = compose_weights(pw₂, m_)
    pₙ = compose_weights(pw₁, n)
    pₙ_ = compose_weights(pw₂, n_)
    return pₙ*pₙ_*ρ^s/(pₘ*pₘ_)
end

function construct_singularity_matrix(μ₁::AbstractInvariantMeasure,
                                    μ₂::AbstractInvariantMeasure,
                                    s::Number;
                                    use_strategy_two::Bool = true
                                    )

    # initialise stuff
    S = [([0],[0])] # needs to be a collection of pairs of indices
    f = [false] # S hasn't been processed yet.
    R = Tuple{Vector{Int64}, Vector{Int64}}[] # blank version of S
    @assert μ₁.supp == μ₂.supp "support of measures must match"
    Γ = μ₁.supp
    pw₁ = μ₁.weights
    G₁ = μ₁.symmetries
    pw₂ = μ₂.weights
    G₂ = μ₂.symmetries
    M = length(Γ.ifs)
    A = zeros(1,1)
    B = zeros(1,1)
    L = zeros(1) # constant vector of log terms, only non-zero when s=0

    pw₁ == pw₂ ? fubuni_flag = true : fubuni_flag = false

    A_rows = 0
    A_cols = 0
    B_rows = 0
    B_cols = 0

    while sum(f)<length(f) # while f contains zeros
        for ∫∫_count = 1:length(S)
            if ~f[∫∫_count]
                a_row = zeros(length(S))
                a_row[∫∫_count] = 1.0
                b_row = zeros(length(R))
                ∫∫_indices = S[∫∫_count]
                if use_strategy_two # subdivide the largest subfractal
                    ∫∫_indices[1]==[0] ? diam_m = Γ.diam : diam_m = Γ.diam*prod([Γ.ifs[m].ρ for m ∈ ∫∫_indices[1]])
                    ∫∫_indices[2]==[0] ? diam_m_ = Γ.diam : diam_m_ = Γ.diam*prod([Γ.ifs[m].ρ for m ∈ ∫∫_indices[2]])
                    if diam_m ≈ diam_m_
                        mrange = 1:M
                        m_range = 1:M
                    elseif diam_m > diam_m_
                        mrange = 1:M
                        m_range = [Int64[]]
                    else
                        mrange = [Int64[]]
                        m_range = 1:M
                    end
                else # subdivide both
                    mrange = 1:M
                    m_range = 1:M
                end
                for m=mrange, m_=m_range
                    mcat = vcat_(∫∫_indices[1],m)
                    mcat_ = vcat_(∫∫_indices[2],m_)
                    is_S_similar, ρ, similar_index = check_for_similar_integrals(Γ, S, mcat, mcat_, G₁, G₂, fubuni_flag)
                    is_ℓ_singular = (check_for_ℓ_singular_integrals(Γ, mcat, mcat_) == 1)

                    # only need to check for R similarities if regular integral.
                    is_singular = is_S_similar || is_ℓ_singular
                    if !is_singular
                        is_R_similar, ρ, similar_index = check_for_similar_integrals(Γ, R, mcat, mcat_, G₁, G₂, fubuni_flag)
                    else
                        is_R_similar = false
                    end
                    is_similar = is_S_similar || is_R_similar
                    if is_similar
                        is_S_similar ? similar_indices = S[similar_index] : similar_indices = R[similar_index]
                        scale_adjust = similar_scaler(ρ, s, similar_indices[1], similar_indices[2], mcat, mcat_, pw₁, pw₂)
                    end

                    if is_ℓ_singular && !is_S_similar # new singularity type
                        push!(S,(mcat,mcat_)) # new type of singularitiy
                        push!(f,false)
                        push!(a_row, -1.0) # increase row by one
                        push!(L,0.0)
                        if A_rows > 0
                            A = hcat(A, zeros(A_rows))
                        end
                    elseif is_S_similar # singular, but seen similar
                        a_row[similar_index] -= scale_adjust
                        if s == 0
                             L[∫∫_count] += μ₁.suppmeasure*μ₂.suppmeasure*log(1/ρ)*prod(pw₁[mcat])*prod(pw₂[mcat_]) # log constant adjustment
                        end
                    elseif is_R_similar # smooth, but seen similar
                        b_row[similar_index] += scale_adjust
                        if s == 0
                            L[∫∫_count] += μ₁.suppmeasure*μ₂.suppmeasure*log(1/ρ)*prod(pw₁[mcat])*prod(pw₂[mcat_]) # log constant adjustment
                        end
                    else # smooth, nothing similar
                        push!(R,(mcat,mcat_))
                        push!(b_row, 1.0) # increase row by one
                        if B_rows > 0
                            B = hcat(B, zeros(B_rows))
                        end
                    end
                end
                # add new row to each matrix
                if A_rows == 0
                    A = reshape(a_row, 1, length(a_row))
                else
                    A = vcat(A, reshape(a_row, 1, length(a_row)))
                end
                if B_rows == 0
                    B = reshape(b_row, 1, length(b_row))
                else
                    B = vcat(B, reshape(b_row, 1, length(b_row)))
                end
                f[∫∫_count] = true
            end
            # update matrix sizes
            A_rows, A_cols = size(A)
            B_rows, B_cols = size(B)
        end
    end
    return A,B,S,R,L
end

"""
    s_energy(Γ::SelfSimilarFractal, s::Number, quad_rule::Function; p₂::Vector{Float64} = getweights(Γ),
     G::Vector{AutomorphicMap}=TrivialGroup(Γ.spatial_dimension),
     G₁::Vector{AutomorphicMap}=TrivialGroup(Γ.spatial_dimension), G₂::Vector{AutomorphicMap}=TrivialGroup(Γ.spatial_dimension))


s is the value in |x-y|⁻ˢ, unless s==0, in which case log|x-y| is used.
p₂ is an optional set of (probability) weights describing an invariant measure of the outer integral.
G₁ and G₂ are groups describing the symmetries of the inner and outer measures respectively.
If G is defined, both measures are assigned this symmetry.
Computes the s-energy of a fractal Γ, using the function quad_rule. This must be of the form:

    quad_rule = (e,j,f) -> I ≈ ∫ₑ∫ⱼ f(x,y) p₁(x)p₂(y)

where A and B are SelfSimilarFractal.
If quad_rule is replaced by some h::Number, the barycentre rule is used with meshwidth h.
"""
function s_energy(  μ₁::AbstractInvariantMeasure,
                    μ₂::AbstractInvariantMeasure,
                    s::Number;
                    h_quad::Float64 = 0.0,
                    N_quad::Int64 = 0,
                    quadrule::Tuple{AbstractArray, AbstractArray, AbstractArray} =
                        getdefault_quad_premap(μ₁, μ₂, h_quad = h_quad, N_quad = N_quad),
                    use_strategy_two::Bool = true
                    )

    @assert μ₁.supp == μ₂.supp "Supports of measures must match."

    # execute main algorithm to express singular integral in terms of smooth integrals
    A,B,_,R,L = construct_singularity_matrix(μ₁, μ₂, s, use_strategy_two = use_strategy_two)
    
    r = zeros(length(R))
    for n in eachindex(r)
        (m,m_) = R[n]
        # map quadrature rule to (m,m') subcomponents
        x, y, w = mapquadrule(μ₁, μ₂, m, m_, quadrule[1], quadrule[2], quadrule[3])
        r[n] = dot(transpose(w), energykernel(s, x, y))
    end
    x = A\(B*r+L)

    return x[1]
end

# common case of using the same measure twice
s_energy(μ, s; vargs...) = s_energy(μ, μ, s; vargs...)

# default to Hausdorff measure when attractor is passed as first argument
s_energy(Γ::AbstractAttractor, s; vargs...) =
    s_energy(HausdorffMeasure(Γ), HausdorffMeasure(Γ), s; vargs...)