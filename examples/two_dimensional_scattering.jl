### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 38d67822-1d98-11ef-1aea-915089b3d1a9
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
end

# ╔═╡ 1115792d-666a-4ca6-9a0a-1a9c6df6b521
using FractalIntegrals, Plots

# ╔═╡ 66ab830e-b2f8-479a-b448-691d549d6f13
import FractalIntegrals: singlelayer_operator_helmholtz,
                            discretise,
                            farfield_pattern,
                            singlelayer_potential_helmholtz

# ╔═╡ becdfa25-8aae-4a4e-a50f-e90d2216a19e
md"""Model incident plane wave $u^i(x) = \mathrm{e}^{\mathrm{i}k d\cdot x}$"""

# ╔═╡ 57f77d24-72e7-41ea-8981-bf4b31c4e1d5
begin
	k = 15
	d = [1, 1] / sqrt(2)
	uⁱ(x) = exp(im*k*((d[1]*x[1])+(d[2]*x[2])))
end

# ╔═╡ 3cd69110-5507-4aaf-9e6b-137e68bfe31e
md"""scattering by fractal $\Gamma$."""

# ╔═╡ 63fb2b9b-1fbb-46d1-bf6f-4dafbd263871
# Γ = getfractal("cantor set", ρ = 0.5);
Γ = getfractal("cantordust", ρ = 1/3);
# Γ = getfractal("koch curve");

# ╔═╡ a0b05573-2cec-4f6b-8842-63613de8a0f8


# ╔═╡ fc699c94-b947-4c3a-8b8c-2b3e75fa7b45
md"""## The Boundary Integral Equation
Define single layer operator on Cantor set $\Gamma$, w.r.t Hausdorff measure:
```math
S_k\phi(x) = \int_\Gamma \Phi(x,y)~\phi(y)~\mathrm{d}\mathcal{H}^d(y),\quad x\in\Gamma.
```
"""

# ╔═╡ 85d36bae-021e-4a52-b767-d65e1f3fc981
Sₖ = singlelayer_operator_helmholtz(Γ, k, ambient_dimension = 2);

# ╔═╡ 5ece5638-417f-4031-a640-a753c307b414
md"""Define RHS data"""

# ╔═╡ 1640b341-3d3f-4a5d-b2c0-d62114ffbf4b
Γ.n == 1 ? f(x) = exp(im*k*d[1]*x) : f(x) = uⁱ(x)

# ╔═╡ 031bfe28-6d37-4e6d-8d27-82ae130e1bff
md"""Consider the BIE
```math
S_k v = f,\quad \text{on }\Gamma.
```
We can solve an approximate version of this for some $v_h\in V_h$, where $V_h$ is a space of piecewise constant basis functions defined on $\Gamma$ with meshwidth $h$. The Galerkin problem is:
```math
(S_k v_h, \varphi_h) = (f, \varphi_h),\quad \varphi\in V_h.
```
"""

# ╔═╡ 76924a03-6acb-4e39-a826-078e684cd9da
Sₖₕ = discretise(Sₖ, h_mesh = 0.01);

# ╔═╡ 2cb80c0b-e9db-4163-8489-6fdbcfca50a3
vₕ = Sₖₕ \ f;

# ╔═╡ 994efbba-b937-41d9-b699-75172008f4ca
md"""## Plotting total field"""

# ╔═╡ c0abac8c-39c2-4a59-bb2d-422555ea0b4f
md"""
Define single layer potential operator on Cantor set $\Gamma$, w.r.t Hausdorff measure:
```math
\mathcal{S}_k\phi(x) = \int_\Gamma \Phi(x,y)~\phi(y)~\mathrm{d}\mathcal{H}^d(y),\quad x\in\mathbb{R}^2.
```
"""

# ╔═╡ d10d1430-92a6-4f6f-a7fe-9f197d0e8d36
slp = singlelayer_potential_helmholtz(vₕ, k, ambient_dimension = 2);

# ╔═╡ f4855596-af3a-4f9c-b1e3-c4a21ba7de9c
begin
δ = 0.01
space_ratio = 0.7  # should be bigger than 0.5

# get centre of fractal
Γ_centre = FractalIntegrals.get_boundingball_centre(Γ)
length(Γ_centre) == 1 ? Γ_centre = [Γ_centre[1],0] : Nothing
x_range = (Γ_centre[1]-space_ratio*Γ.diam):δ:(Γ_centre[1]+space_ratio*Γ.diam)
y_range = (Γ_centre[2]-space_ratio*Γ.diam):δ:(Γ_centre[2]+space_ratio*Γ.diam)
plot_pts = [[x,y] for x in x_range, y in y_range]
end;

# ╔═╡ 7bd71f5c-9fc4-449c-a3eb-08efe3866c40
begin
# create the total field function
uˢ(x) = -slp(x)
u(x) = uⁱ(x) + uˢ(x)

# plot the total field
heatmap(x_range,y_range,transpose(real.(u.(plot_pts))), aspect_ratio = 1, colormap=:jet)
# plot the fractal on top
plot!(Γ,color=:black, markersize = 0.3)
end

# ╔═╡ 4633dc4f-7451-46a4-863c-f8dd2cc40d86
md"""## Plotting far field pattern"""

# ╔═╡ bd4de407-a217-4d6a-99e0-b98f35efca37
θ = 0:0.01:2π

# ╔═╡ c6ffa209-db0f-4ff5-84e9-84ad165fc3db
ffp = farfield_pattern(ϕ, k)

# ╔═╡ 8f26187a-8dbe-4e3d-8b9b-322d2296a8da
plot(θ, [real.(ffp.(θ)), imag.(ffp.(θ))], labels = ["realffp" "imagffp"])

# ╔═╡ Cell order:
# ╠═38d67822-1d98-11ef-1aea-915089b3d1a9
# ╠═1115792d-666a-4ca6-9a0a-1a9c6df6b521
# ╠═66ab830e-b2f8-479a-b448-691d549d6f13
# ╟─becdfa25-8aae-4a4e-a50f-e90d2216a19e
# ╠═57f77d24-72e7-41ea-8981-bf4b31c4e1d5
# ╟─3cd69110-5507-4aaf-9e6b-137e68bfe31e
# ╠═63fb2b9b-1fbb-46d1-bf6f-4dafbd263871
# ╠═a0b05573-2cec-4f6b-8842-63613de8a0f8
# ╟─fc699c94-b947-4c3a-8b8c-2b3e75fa7b45
# ╠═85d36bae-021e-4a52-b767-d65e1f3fc981
# ╟─5ece5638-417f-4031-a640-a753c307b414
# ╠═1640b341-3d3f-4a5d-b2c0-d62114ffbf4b
# ╟─031bfe28-6d37-4e6d-8d27-82ae130e1bff
# ╠═76924a03-6acb-4e39-a826-078e684cd9da
# ╠═2cb80c0b-e9db-4163-8489-6fdbcfca50a3
# ╟─994efbba-b937-41d9-b699-75172008f4ca
# ╟─c0abac8c-39c2-4a59-bb2d-422555ea0b4f
# ╠═d10d1430-92a6-4f6f-a7fe-9f197d0e8d36
# ╠═f4855596-af3a-4f9c-b1e3-c4a21ba7de9c
# ╠═7bd71f5c-9fc4-449c-a3eb-08efe3866c40
# ╟─4633dc4f-7451-46a4-863c-f8dd2cc40d86
# ╠═bd4de407-a217-4d6a-99e0-b98f35efca37
# ╠═c6ffa209-db0f-4ff5-84e9-84ad165fc3db
# ╠═8f26187a-8dbe-4e3d-8b9b-322d2296a8da
