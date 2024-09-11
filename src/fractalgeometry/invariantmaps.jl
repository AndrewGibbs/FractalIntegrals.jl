InvariantMap(δ, A) = Similarity(1.0, δ, A)

# --------------------------- rotation and reflection matrices ----------------------------- #
rotationmatrix2d(θ::T) where T = SMatrix{2,2,T}([cos(θ) -sin(θ); sin(θ) cos(θ)])
reflectionmatrix2d(θ::T) where T = SMatrix{2,2,T}([cos(2θ) sin(2θ); sin(2θ) -cos(2θ)])

# --------------------------- Symmetric groups --------------------------------------------#
function get_group_operations2d(num_rotations::Integer,
                                reflections::AbstractVector{T},
                                centre::AbstractVector{T}) where T<:Number
    δθ = T(2π) / num_rotations
    num_reflections = length(reflections)
    
    ivmaps = Vector{Similarity{2,T}}(undef, num_rotations + num_reflections)

    counter = 0
    for θ in 0:δθ:(2π-δθ)
        counter += 1
        rotmat = rotationmatrix2d(θ)
        ivmaps[counter] = InvariantMap(centre - rotmat * centre, rotmat)
    end
    for ϑ in reflections
        counter += 1
        refmat = reflectionmatrix2d(ϑ)
        ivmaps[counter] = InvariantMap(centre - refmat * centre, refmat)
    end
    # the union of rotations and reflections gives all compositions
    return ivmaps
end

function DihedralGroup( T::Type,
                        n::Integer;
                        centre::AbstractVector{<:Number} = zeros(T, 2),
                        angle_offest::Number = zero(T))
    if mod(n,2) == 0 # for even n, want to double frequency of reflections, but restrict to [0,π]
        δθ = T(π)/n
        reflections = (0:δθ:(T(π)-δθ)) .+ angle_offest
    else
        δθ = 2π/n # for odd n, split [0,2π] into n segments
        reflections = (0:δθ:(2T(π)-δθ)) .+ angle_offest
    end
    return get_group_operations2d(n, reflections, centre)
end

DihedralGroup(n; varags...) = DihedralGroup(Float64, n; varags...)

function trivialgroup(  T::Type,
                        n::Integer)
    if n==1
        tg = Similarity(one(T), zero(T))
    else
        tg = Similarity(one(T), zeros(T, n))
        # was 'TranslatingSimilarity', but this type has been removed
    end
    return [tg]
end
trivialgroup(n::Integer) = trivialgroup(Float64, n)