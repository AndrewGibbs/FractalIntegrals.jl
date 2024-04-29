# needs to be lower case!
function singlelayer_operator_helmholtz(  μ::AbstractInvariantMeasure,
                                        k::Number;
                                        ambient_dimension::Integer = μ.supp.n
    )
    if ambient_dimension == 2     
        K = SelfAdjointSingularIntegralOperator(μ, #fractal domain
        (x,y)->helmholtzkernel2d(k,x,y), # Hankel function
        (x,y)->helmholtzkernel2d_lipschitzpart(k,x,y), # kernel minus singularity
        0.0, # strength of singularity, corresponding to log singularity
        ComplexF64(-1/(2π)), # scaling of singularity
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = SelfAdjointSingularIntegralOperator(μ, #fractal domain
            (x,y)->helmholtzkernel3d(k,x,y), # Green's function
            (x,y)->helmholtzkernel3d_lipschitzpart(k,x,y), # kernel minus singularity
            1.0, # strength of singularity, corresponding to 1/|x-y|
            ComplexF64(1/(4π)), # scaling of singularity
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