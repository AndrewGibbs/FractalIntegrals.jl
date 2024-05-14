
function singlelayer_operator_laplace(μ::AbstractInvariantMeasure{
                                        <:AbstractAttractor{<:Any, R}
                                        },
                                    ambient_dimension::Integer = μ.supp.n
                                    ) where {R <: Real}
    if ambient_dimension == 2
        K = SingularIntegralOperator(Γ, #fractal domain
        (x,y)->-1/(2π)*Φₜ(0.0,x,y), # log kernel
        (x,y)->zero_kernel(x,y), # kernel minus singularity
        0.0, # strength of singularity, corresponding to log singularity
        -1/(2π), # scaling of singularity
        true, #self-adjoint
        0.0 #wavenumber
        )
    elseif ambient_dimension == 3
        K = SIO{Ω}(Γ, #fractal domain
        (x,y)-> 1/(4π)*Φₜ(1.0,x,y), # Green's function
        (x,y)-> zero_kernel(x,y), # kernel minus singularity
        1.0, # strength of singularity, corresponding to 1/|x-y|
        1/(4π), # scaling of singularity
        true, #self-adjoint
        0.0 #wavenumber
        )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

function singlelayer_operator_helmholtz(μ::AbstractInvariantMeasure{
                                                <:AbstractAttractor{<:Any, R}
                                            },
                                        k::Number;
                                        ambient_dimension::Integer = μ.supp.n
    ) where {R <: Real}
    if ambient_dimension == 2     
        K = OscillatorySingularIntegralOperator(μ, #fractal domain
        (x,y) -> energykernel(0, x, y), # Hankel function
        (x,y) -> zero_kernel(x,y), # kernel minus singularity
        R(0.0), # strength of singularity, corresponding to log singularity
        Complex{R}(-1/(2π)), # scaling of singularity
        true, # self-adjoint
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = OscillatorySingularIntegralOperator(μ, #fractal domain
            (x,y) -> energykernel(1,x,y), # Green's function
            (x,y) -> zero_kernel(x,y), # kernel minus singularity
            R(1.0), # strength of singularity, corresponding to 1/|x-y|
            ComplexF64{R}(1/(4π)), # scaling of singularity
            true, # self-adjoint
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