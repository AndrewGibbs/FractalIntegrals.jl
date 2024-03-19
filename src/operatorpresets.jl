function SingleLayerOperatorHelmholtz(  μ::AbstractInvariantMeasure,
                                        k::Number;
                                        ambient_dimension::Integer = μ.supp.n
    )
    if ambient_dimension == 2     
        K = SeparableIntegralOperator(μ, #fractal domain
        (x,y)->HelhmoltzGreen2D(k,x,y), # Hankel function
        (x,y)->HelhmoltzGreen2D_Lipschitz_part(k,x,y), # kernel minus singularity
        0.0, # strength of singularity, corresponding to log singularity
        -1/(2π), # scaling of singularity
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = SeparableIntegralOperator(μ, #fractal domain
            (x,y)->HelhmoltzGreen3D(k,x,y), # Green's function
            (x,y)->HelhmoltzGreen3D_Lipschitz_part(k,x,y), # kernel minus singularity
            1.0, # strength of singularity, corresponding to 1/|x-y|
            1/(4π), # scaling of singularity
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end