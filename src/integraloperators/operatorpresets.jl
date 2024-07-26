function check_ambient_dimension(ambient_dimension)
    if ambient_dimension == 1
        @warn "Cannot have ambient_dimension=1, assuming two-dimensional 'screen' problem"
        ambient_dimension = 2
    end
    return ambient_dimension
end

function singlelayer_operator_laplace(μ::AbstractInvariantMeasure{
                                        <:AbstractAttractor{R, <:Any}
                                        };
                                    ambient_dimension::Integer = μ.supp.n
                                    ) where {R <: Real}

    ambient_dimension = check_ambient_dimension(ambient_dimension)

    if ambient_dimension == 2     
        K = SingularIntegralOperator(μ, #fractal domain
        (x,y) -> energykernel(0, x, y), # Hankel function
        (x,y) -> zero_kernel(x,y), # kernel minus singularity
        R(0.0), # strength of singularity, corresponding to log singularity
        R(-1/(2π)), # scaling of singularity
        true, # self-adjoint
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = SingularIntegralOperator(μ, #fractal domain
            (x,y) -> energykernel(1,x,y), # Green's function
            (x,y) -> zero_kernel(x,y), # kernel minus singularity
            R(1.0), # strength of singularity, corresponding to 1/|x-y|
            R(1/(4π)), # scaling of singularity
            true, # self-adjoint
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end


# Hausdorff default
singlelayer_operator_laplace(Γ::AbstractAttractor; vargs...) = 
    singlelayer_operator_laplace(HausdorffMeasure(Γ); vargs...)


function singlelayer_operator_helmholtz(μ::AbstractInvariantMeasure{
                                                <:AbstractAttractor{R, <:Any}
                                            },
                                        k::Number;
                                        ambient_dimension::Integer = μ.supp.n
    ) where {R <: Real}

    ambient_dimension = check_ambient_dimension(ambient_dimension)

    if ambient_dimension == 2     
        K = OscillatorySingularIntegralOperator(μ, #fractal domain
        (x,y) -> helmholtzkernel2d(k, x, y), # Hankel function
        (x,y) -> helmholtzkernel2d_lipschitzpart(k,x,y), # kernel minus singularity
        R(0.0), # strength of singularity, corresponding to log singularity
        Complex{R}(-1/(2π)), # scaling of singularity
        true, # self-adjoint
        k
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = OscillatorySingularIntegralOperator(μ, #fractal domain
            (x,y) -> helmholtzkernel3d(k,x,y), # Green's function
            (x,y) -> helmholtzkernel3d_lipschitzpart(k,x,y), # kernel minus singularity
            R(1.0), # strength of singularity, corresponding to 1/|x-y|
            ComplexF64{R}(1/(4π)), # scaling of singularity
            true, # self-adjoint
            k
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

# Hausdorff default
singlelayer_operator_helmholtz(Γ::AbstractAttractor, k; vargs...) = 
    singlelayer_operator_helmholtz(HausdorffMeasure(Γ), k; vargs...)

# @default_to_hausdorff singlelayer_operator_laplace
# @default_to_hausdorff singlelayer_operator_helmholtz

# for preset_fn ∈ [SingleLayerOperatorHelmholtz]
#     preset_fn(Γ::Attractor, args...; vargs...) = preset_fn(HausdorffMeasure(Γ), args...; vargs...)
# end