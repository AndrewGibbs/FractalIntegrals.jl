# subdivide attractor and get array of VectorIndex 
function subdivide_indices( Γ::AbstractAttractor{N, M, T},
                            h::Real;
                            max_num_indices = Inf
                            ) where {N, M, T}


    @assert (h>0 || max_num_indices<Inf
    ) "either meshwidth must be positive, or max_num_indices must be finite"

    Lₕ = [zero(VectorIndex{M, Int64})]#[VectorIndex{M}([0])]
    ρs = [s.ρ for s in Γ.ifs]
    diams = [diam(Γ)]

    keep_subdividing = true
    while keep_subdividing && length(Lₕ)<max_num_indices
        split_vecs = Int64[]
        keep_subdividing = false
        for j in eachindex(Lₕ)
            if diams[j] ≥ h
                keep_subdividing = true
                # append new vector index onto end of index set
                append!(Lₕ, split(Lₕ[j]))
                # similar for diameter logging
                append!(diams, diams[j].*ρs)

                # log current indices to delete
                push!(split_vecs,j)
            end
        end

        # delete vectors and radii from split components
        deleteat!(Lₕ, split_vecs)
        deleteat!(diams, split_vecs)
    end

    return Lₕ
end

function grade_mesh( Γ::AbstractAttractor{N, M, T},
                        subdivide_if_true_fn::Function;
                        max_num_indices = 1e6
                        ) where {N, M, T}

    Lₕ = [zero(VectorIndex{M, Int64})]
    mesh = [Γ]

    keep_subdividing = true
    while keep_subdividing && (length(Lₕ) < max_num_indices)
        keep_subdividing = false
        for (j, Γₘ) in enumerate(mesh)
            if subdivide_if_true_fn(Γₘ)
                # add new index vectors and mesh elements
                append!(Lₕ, split(Lₕ[j]))
                append!(mesh, split(mesh[j]))

                # delete vectors and mesh els from split components
                deleteat!(Lₕ, j)
                deleteat!(mesh, j)

                # break outer for loop and restart for loop (via while loop) 
                # ... to avoid messing with indices after changing vector size
                keep_subdividing = true
                break
            end
        end
    end

    return mesh, Lₕ
end

grade_towards_points_fn(Γ::AbstractAttractor, x, C::Number) = dist⁻(Γ, x)/diam(Γ) < C
grade_mesh_towards_point(Γ::AbstractAttractor, x; C::Number=1) = 
    grade_mesh(Γ, γ -> grade_towards_points_fn(γ, x, C))