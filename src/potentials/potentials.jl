struct Potential{P<:Projection,
                F<:Function,
                Q<:AbstractArray{<:QuadStruct}}
    density::P # measure, support etc contained in here
    kernel::F
    quadrules::Q
end

function (pot::Potential)(x)
    val = zero(typeof(x))
    @inbounds @simd for n in 1:length(pot.density.basis)
        @fastmath val += (pot.density.coeffs[n] * (
                transpose(pot.quadrules[n].weights) *
                (pot.kernel.(x, pot.quadrules[n].nodes) .* pot.density.basis[n].(pot.quadrules[n].nodes))
                ))
    end
    return val
end

# function parapot(pot::Potential, x::AbstractArray{T}) where T<:Number
#     tasks = map(1:length(pot.density.basis)) do n
#         @spawn pot.density.coeffs[n] * 
#             dot(pot.quadrules[n][2],
#             pot.kernel(x, pot.quadrules[n][1]) .* pot.density.basis[n](pot.quadrules[n][1]))
#     end
#     return sum(fetch.(tasks))
# end

# function (pot::Potential)(X::AbstractArray{<:AbstractArray{<:Number}})
#     U = Matrix{ComplexF64}(undef,size(X))
#     # tasks = chunks(X, n=10*nthreads())
#     tasks = collect(Iterators.partition(eachindex(X), div(length(X), 10*nthreads())))
#     # Xtasks = [@views X[ts] for ts in tasks]
#     @sync for n in 1:length(tasks)
#         @spawn begin
#             Xloc = X[tasks[n]]
#             @inbounds U[tasks[n]] = pot.(Xloc)
#         end
#     end
#     # Threads.@threads for j in eachindex(X)
#     #     U[j] = pot(X[j])
#     # end
#     return U
# end

include("potentialpresets.jl")