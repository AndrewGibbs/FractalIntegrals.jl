# needs to be lower case!
function singlelayer_operator_helmholtz(μ::AbstractInvariantMeasure{
                                                <:AbstractAttractor{<:Any, R}
                                            },
                                        k::Number;
                                        ambient_dimension::Integer = μ.supp.n
    ) where {R <: Real}
    if ambient_dimension == 2     
        K = OscillatorySingularIntegralOperator(μ, #fractal domain
        (x,y)->helmholtzkernel2d(R(k),x,y), # Hankel function
        (x,y)->helmholtzkernel2d_lipschitzpart(R(k),x,y), # kernel minus singularity
        R(0.0), # strength of singularity, corresponding to log singularity
        Complex{R}(-1/(2π)), # scaling of singularity
        true, # self-adjoint
        R(k) # store wavenumber, useful for automating stuff later
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = OscillatorySingularIntegralOperator(μ, #fractal domain
            (x,y)->helmholtzkernel3d(R(k),x,y), # Green's function
            (x,y)->helmholtzkernel3d_lipschitzpart(R(k),x,y), # kernel minus singularity
            R(1.0), # strength of singularity, corresponding to 1/|x-y|
            ComplexF64{R}(1/(4π)), # scaling of singularity
            true, # self-adjoint
            R(k), # store wavenumber, useful for automating stuff later
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

singlelayer_operator_helmholtz(Γ::AbstractAttractor, args...; vargs...) = 
    singlelayer_operator_helmholtz(HausdorffMeasure(Γ), args...; vargs...)

# for preset_fn ∈ [SingleLayerOperatorHelmholtz]
#     preset_fn(Γ::Attractor, args...; vargs...) = preset_fn(HausdorffMeasure(Γ), args...; vargs...)
# end