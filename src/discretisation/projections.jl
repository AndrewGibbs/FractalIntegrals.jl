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

function Base.:\(op::FractalOperator, f::Function)
    discop = discretise(op)
    fₕ = project(   discop.ip, # inner product
                    discop.basis, # basis
                    f)
    coeffs = discop.galerkinmatrix \ fₕ.coeffs
    return Projection(discop.basis, coeffs)
end