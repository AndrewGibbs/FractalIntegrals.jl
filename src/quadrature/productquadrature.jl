# the below function should be in a more generic location.
# also, there should be a fast version for scalar w, i.e. when all weights are equal
function combine_quadrules(x1, w1, x2, w2)
    return repeat(x1; inner = length(w2)),
    repeat(x2; outer = length(w1)),
    repeat(w1; inner = length(w2)) .* repeat(w2; outer = length(w1))
end

function combine_quadrules(
    x1::AbstractArray{<:Union{SVector,Number}},
    w1::AbstractArray{<:Number},
    x2::AbstractArray{<:Union{SVector,Number}},
    w2::AbstractArray{<:Number},
)
    N1 = length(w1)
    N2 = length(w2)
    N = N1 * N2
    X1 = Vector{eltype(x1)}(undef, N)
    X2 = Vector{eltype(x2)}(undef, N)
    W1 = Vector{eltype(w1)}(undef, N)
    W2 = Vector{eltype(w2)}(undef, N)

    @simd for n2 = 1:N2
        X1[((n2-1)*N1+1):(n2*N1)] .= x1
        W1[((n2-1)*N1+1):(n2*N1)] .= w1
    end

    @simd for n2 = 1:N2
        @inbounds for n2_ = ((n2-1)*N1+1):(n2*N1)
            X2[n2_] = x2[n2]
        end
        W2[((n2-1)*N1+1):(n2*N1)] .= w2[n2]
    end

    return X1, X2, W1 .* W2
end

function combine_quadrules(
    x1::AbstractArray{<:Union{SVector,Number}},
    w1::Number,
    x2::AbstractArray{<:Union{SVector,Number}},
    w2::Number,
)
    N1 = length(x1)
    N2 = length(x2)
    N = N1 * N2
    X1 = Vector{eltype(x1)}(undef, N)
    X2 = Vector{eltype(x2)}(undef, N)

    @simd for n1 = 1:N1
        X1[((n1-1)*N2+1):(n1*N2)] .= x1
    end

    @simd for n2 = 1:N2
        @inbounds for n2_ = ((n2-1)*N1+1):(n2*N1)
            X2[n2_] = x2[n2]
        end
    end

    return X1, X2, fill(w1 * w2, N)
end

function combine_quadrules(q1::QuadStruct, q2::QuadStruct)
    return combine_quadrules(q1.nodes, q1.weights, q2.nodes, q2.weights)
end

function tensor_product(t...)
    # get necessary dimension sizes
    num_output_vecs = length(t)
    final_length = prod(length.([t...]))
    # apply tensor product, and flatten
    t_long = reshape(collect(Iterators.product(t...)), final_length, 1)
    # restructure and return
    return ([t_long[j][n] for j = 1:final_length] for n = 1:num_output_vecs)
end

function combine_quadrules(queple::QuadStruct...)# where N
    nodes_tp = tensor_product((q.nodes for q in queple)...)
    weights_tp = tensor_product((q.weights for q in queple)...)
    return nodes_tp..., reduce(.*, weights_tp)
end

function tensor_product_quadrature(quad_fn::Function)
    return function (inputs...)
        μple = inputs[1:(end-1)]
        param = inputs[end]
        return (combine_quadrules((QuadStruct(quad_fn(μ, param)...) for μ in μple)...))
    end
end

macro tensorquad(quadrule)
    new_name = Symbol("tensor_", quadrule)
    quote
        $(esc(new_name)) = tensor_product_quadrature($(esc(quadrule)))
    end
end
