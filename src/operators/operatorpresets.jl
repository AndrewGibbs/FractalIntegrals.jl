function check_ambient_dimension(ambient_dimension)
    if ambient_dimension == 1
        @warn "Cannot have ambient_dimension=1, assuming two-dimensional 'screen' problem"
        ambient_dimension = 2
    end
    return ambient_dimension
end

function singlelayer_operator_laplace(μ::AbstractInvariantMeasure{
                                        N, <:Any, R, <:AbstractAttractor
                                        };
                                    ambient_dimension::Integer = N
                                    ) where {N, R <: Real}

    ambient_dimension = check_ambient_dimension(ambient_dimension)

    if ambient_dimension == 2     
        K = SeparableIntegralOperator(μ, #fractal domain
        (x, y) -> energykernel(0, x, y), # Hankel function
        (x, y) -> zero_kernel(x, y), # kernel minus singularity
        zero(R), # strength of singularity, corresponding to log singularity
        R(-1/(2π)), # scaling of singularity
        true, # self-adjoint
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = SeparableIntegralOperator(μ, #fractal domain
            (x, y) -> energykernel(1, x, y), # Green's function
            (x, y) -> zero_kernel(x, y), # kernel minus singularity
            one(R), # strength of singularity, corresponding to 1/|x-y|
            R(1/(4π)), # scaling of singularity
            true, # self-adjoint
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

# Hausdorff default
@hausdorffdefault singlelayer_operator_laplace

function singlelayer_operator_helmholtz(μ::AbstractInvariantMeasure{
                                            N, <:Any, R, <:AbstractAttractor
                                            },
                                        k::Number;
                                        ambient_dimension::Integer = N
    ) where {N, R <: Real}

    ambient_dimension = check_ambient_dimension(ambient_dimension)

    if ambient_dimension == 2     
        K = OscillatorySeparableIntegralOperator(μ, #fractal domain
        (x,y) -> helmholtzkernel2d(k, x, y), # Hankel function
        (x,y) -> helmholtzkernel2d_lipschitzpart(k, x, y), # kernel minus singularity
        zero(R), # strength of singularity, corresponding to log singularity
        Complex{R}(-1/(2π)), # scaling of singularity
        true, # self-adjoint
        R(k)
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = OscillatorySeparableIntegralOperator(μ, #fractal domain
            (x,y) -> helmholtzkernel3d(k,x,y), # Green's function
            (x,y) -> helmholtzkernel3d_lipschitzpart(k,x,y), # kernel minus singularity
            one(R), # strength of singularity, corresponding to 1/|x-y|
            ComplexF64{R}(1/(4π)), # scaling of singularity
            true, # self-adjoint
            R(k)
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

# Hausdorff default
@hausdorffdefault singlelayer_operator_helmholtz