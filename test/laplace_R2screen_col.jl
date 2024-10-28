import FractalIntegrals: get_barycentre, barycentre_quadrule, diam

# anti-derivative of Laplace kernel
function log_antideriv(x::Real)
    if x>0
        return -x*(log(x)-1)/(2π) # integrate by parts (I think)
    elseif x==0
        return 0.0 # take limit as x→0
    else
        error("input to log anti-deriv must be positive")
    end
end

# exact value of integral of ∫_[a,b] log(|x-y|)dy
function laplace_collocation_integral(a::Real, b::Real, x::Real)
    if b ≤ x
        return log_antideriv(x-a) - log_antideriv(x-b)
    elseif x ≤ a
        return log_antideriv(b-x) - log_antideriv(a-x)
    else #a≤x≤b
        return laplace_collocation_integral(a, x, x) + laplace_collocation_integral(x, b, x)
    end
end

# returns the exact (up to machine accuracy) collocation matrix for this special case
function get_exact_laplace_screen_col_matrix(basis, collocation_points)
    col_mat_exact = zeros((length(collocation_points),length(basis)))
    for (n,ϕ) in enumerate(basis)
        for (m,x) in enumerate(collocation_points)
            # get integration interval [a, b]
            a = get_barycentre(ϕ.measure) - diam(ϕ.measure)/2
            b = get_barycentre(ϕ.measure) + diam(ϕ.measure)/2
            col_mat_exact[m,n] = laplace_collocation_integral(a,b,x.node)
        end
    end
    return col_mat_exact
end