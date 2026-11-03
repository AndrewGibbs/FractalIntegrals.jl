function check_ambient_dimension(ambient_dimension)
    if ambient_dimension == 1
        @warn "Cannot have ambient_dimension=1, assuming two-dimensional 'screen' problem"
        ambient_dimension = 2
    end
    return ambient_dimension
end

## --------------------- single layer opeartor Helmholtz -------------------- ##

"""
    singlelayer_operator_laplace(μ::AbstractInvariantMeasure; ambient_dimension)
    singlelayer_operator_laplace(Γ::AbstractAttractor; ambient_dimension)
    singlelayer_operator_laplace(μple::Tuple; ambient_dimension)

Construct the Laplace single-layer integral operator on a fractal measure support
``\\Gamma``, defined by

```math
\\mathcal{K}\\,\\phi\\,(x) = \\int_\\Gamma \\Phi(x,y)\\,\\phi(y)\\,\\mathrm{d}\\mu(y),
\\quad x \\in \\Gamma,
```

where ``\\Phi`` is the free-space Green's function of the Laplacian:

- **2D**: ``\\Phi(x,y) = -\\tfrac{1}{2\\pi}\\log|x-y|``
- **3D**: ``\\Phi(x,y) = \\tfrac{1}{4\\pi}|x-y|^{-1}``

The kernel is decomposed into singular and smooth parts as a
`SeparableIntegralOperator`, enabling efficient quadrature on fractal domains.

When a `Tuple` of measures (or attractors) is passed, a `BlockOperator` is
returned representing the operator on the union of supports.

# Arguments
- `μ`: An `AbstractInvariantMeasure` supported on ``\\Gamma``.
- `Γ`: An `AbstractAttractor` (uses Hausdorff measure by default).
- `μple`: A `Tuple` of measures for a multi-component domain.
- `ambient_dimension`: The ambient space dimension (default: inferred from `μ`).
  Note: `ambient_dimension=1` is not valid and will be promoted to `2`.

# Examples
```julia
Γ = getfractal("cantorset")
K = singlelayer_operator_laplace(Γ)
```
"""
function singlelayer_operator_laplace(μ::AbstractInvariantMeasure{
                                        N, <:Any, R, <:AbstractAttractor
                                        };
                                    ambient_dimension::Integer = N
                                    ) where {N, R <: Real}

    ambient_dimension = check_ambient_dimension(ambient_dimension)

    if ambient_dimension == 2     
        K = SeparableIntegralOperator(μ, #fractal domain
        (x, y) -> -energykernel(0, x, y)/(2π), # log function
        (x, y) -> zero_kernel(x, y), # kernel minus singularity
        zero(R), # strength of singularity, corresponding to log singularity
        R(-1/(2π)), # scaling of singularity
        true, # self-adjoint
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = SeparableIntegralOperator(μ, #fractal domain
            (x, y) -> energykernel(1, x, y)/(4π), # Green's function
            (x, y) -> zero_kernel(x, y), # kernel minus singularity
            one(R), # strength of singularity, corresponding to 1/|x-y|
            R(1/(4π)), # scaling of singularity
            true, # self-adjoint
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

function dom2codom_singlelayer_operator_laplace(μ₁::M1,
                                                μ₂::M2;
                                                ambient_dimension = max(get_ambient_dimension.(μ₁, μ₂)...),
                                                varargs...) where {M1<:AbstractInvariantMeasure, M2<:AbstractInvariantMeasure}
    T = eltype(eltype(μ₁))
    if μ₁ == μ₂
        singlelayer_operator_laplace(μ₁; varargs...)
    elseif ambient_dimension == 2     
        Φ = (x, y) -> -energykernel(0, x, y)/(2π)
        K = SmoothIntegralOperator{M1, M2, typeof(Φ), T}(
                            μ₁,
                            μ₂, #domain -> codomain
                            Φ, # Hankel function
                            false
                        )
    elseif ambient_dimension == 3
        #3D Helmholtz case  
        Φ = (x, y) -> energykernel(1, x, y)/(4π)
        K = SmoothIntegralOperator{M1, M2, typeof(Φ), T}(
                            μ₁,
                            μ₂, #domain -> codomain
                            Φ, # Hankel function
                            false
                        )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

function singlelayer_operator_laplace(μple::Tuple{Vararg{AbstractInvariantMeasure}};
                                    ambient_dimension = max.(get_ambient_dimension.(μple)...),
                                    varargs...)

    # first make sure all supports have the same ambient dimension
    μple = embed_into_same_dimension(μple)
                        
    # (note that this may be less than the ambient dimension specified above)

    # now define the block of operators
    op_block =
        [dom2codom_singlelayer_operator_laplace(μ₁, μ₂;
                                                ambient_dimension = ambient_dimension,
                                                varargs...)
                                                for μ₂ in μple, μ₁ in μple
        ]
    return BlockOperator(∪(μple...), op_block, true)
end

# Hausdorff defaults
@hausdorffdefault singlelayer_operator_laplace
singlelayer_operator_laplace(Γ_tuple::Tuple{Vararg{AbstractAttractor}}; varargs...) =
    singlelayer_operator_laplace(HausdorffMeasure.(Γ_tuple); varargs...)


## --------------------- single layer opeartor Helmholtz -------------------- ##

"""
    singlelayer_operator_helmholtz(μ::AbstractInvariantMeasure, k; ambient_dimension)
    singlelayer_operator_helmholtz(Γ::AbstractAttractor, k; ambient_dimension)
    singlelayer_operator_helmholtz(μple::Tuple, k; ambient_dimension)

Construct the Helmholtz single-layer integral operator with wavenumber `k` on a
fractal measure support ``\\Gamma``, defined by

```math
\\mathcal{K}\\,\\phi\\,(x) = \\int_\\Gamma \\Phi_k(x,y)\\,\\phi(y)\\,\\mathrm{d}\\mu(y),
\\quad x \\in \\Gamma,
```

where ``\\Phi_k`` is the free-space Green's function of the Helmholtz equation
``(\\Delta + k^2)u = 0``:

- **2D**: ``\\Phi_k(x,y) = \\tfrac{\\mathrm{i}}{4} H_0^{(1)}(k|x-y|)``
- **3D**: ``\\Phi_k(x,y) = \\tfrac{e^{\\mathrm{i}k|x-y|}}{4\\pi|x-y|}``

The oscillatory kernel is decomposed into a log-singular (2D) or ``1/r``-singular
(3D) part and a smooth Lipschitz remainder, stored as an
`OscillatorySeparableIntegralOperator`.

When a `Tuple` of measures (or attractors) is passed, a `BlockOperator` is
returned representing the operator on the union of supports.

# Arguments
- `μ`: An `AbstractInvariantMeasure` supported on ``\\Gamma``.
- `Γ`: An `AbstractAttractor` (uses Hausdorff measure by default).
- `μple`: A `Tuple` of measures for a multi-component domain.
- `k`: The wavenumber ``k > 0``.
- `ambient_dimension`: The ambient space dimension (default: inferred from `μ`).
  Note: `ambient_dimension=1` is not valid and will be promoted to `2`.

# Examples
```julia
Γ = getfractal("cantordust")
k = 10.0
K = singlelayer_operator_helmholtz(Γ, k)
```
"""
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
            Complex{R}(1/(4π)), # scaling of singularity
            true, # self-adjoint
            R(k)
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end


function dom2codom_singlelayer_operator_helmholtz(  μ₁::M1,
                                                    μ₂::M2,
                                                    k::Number;
            ambient_dimension = max(get_ambient_dimension.(μ₁, μ₂)...),
            varargs...) where {M1<:AbstractInvariantMeasure, M2<:AbstractInvariantMeasure}
    T = eltype(eltype(μ₁))
    if μ₁ == μ₂
        singlelayer_operator_helmholtz(μ₁, k; varargs...)
    elseif ambient_dimension == 2     
        Φ = (x, y) -> helmholtzkernel2d(k, x, y)
        K = SmoothIntegralOperator{M1, M2, typeof(Φ), T}(
        μ₁,
        μ₂, #domain -> codomain
        Φ, # Hankel function
        false
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case  
        Φ = (x, y) -> helmholtzkernel3d(k, x, y)
        K = SmoothIntegralOperator{M1, M2, typeof(Φ), T}(
        μ₁,
        μ₂, #domain -> codomain
        Φ, # Hankel function
        false
        )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

function singlelayer_operator_helmholtz(μple::Tuple{Vararg{AbstractInvariantMeasure}},
                                        k::Number;
                                        ambient_dimension = max.(get_ambient_dimension.(μple)...),
                                        varargs...)

    # first make sure all supports have the same ambient dimension
    μple = embed_into_same_dimension(μple)

    # (note that this may be less than the ambient dimension specified above)

    # now define the block of operators
    op_block =
    [dom2codom_singlelayer_operator_helmholtz(μ₁, μ₂, k;
        ambient_dimension = ambient_dimension,
        varargs...)
        for μ₂ in μple, μ₁ in μple
    ]
    return BlockOperator(∪(μple...), op_block, true)
end


# Hausdorff defaults
@hausdorffdefault singlelayer_operator_helmholtz
singlelayer_operator_helmholtz(Γ_tuple::Tuple{Vararg{AbstractAttractor}}, k::Number; varargs...) =
    singlelayer_operator_helmholtz(HausdorffMeasure.(Γ_tuple), k; varargs...)

# create alias for SLO
newton_operator_helmholtz = singlelayer_operator_helmholtz