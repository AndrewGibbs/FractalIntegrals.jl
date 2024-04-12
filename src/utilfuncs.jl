function promote_parameter(v,T)
    u = similar(v,T)
    for j in eachindex(u)
        u[j] = T(v[j])
    end
    return u
end