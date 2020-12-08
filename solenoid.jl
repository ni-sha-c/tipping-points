### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ b49dd7c8-3772-11eb-3414-bb1f9e84b748
begin
	using LaTeXStrings
	using Plots
	using PlutoUI
end

# ╔═╡ 7abc4f04-3772-11eb-1c9b-712985ec2af7
md"""
### Dependencies
"""

# ╔═╡ f8e62384-3772-11eb-23fb-41462d88e988
function solenoid(x, s)
	
	s₀, s₁ = s[1], s[2]
	
	r = sqrt(x[1]*x[1] + x[2]*x[2])
	θ = (atan(x[2],x[1]) + 2π) % (2π) 
	z = x[3]
	
	r = s₀ + (r - s₀)/s₁ + cos(θ)/2
	θ = 2θ
	z = z/s₁ + sin(θ)/2
	
	xn = similar(x)
	xn[1] = r*cos(θ)
	xn[2] = r*sin(θ)
	xn[3] = z
	
	return xn
end

# ╔═╡ 9d23a890-3773-11eb-1d5e-295158261eea
function run(method, x₀, s, n_steps)
	n = 0
	orbit = [x₀]
	x₁ = copy(x₀)
	while n < n_steps
		x₁ = method(x₁, s)
		push!(orbit, x₁)
		n += 1
	end
	return orbit
end

# ╔═╡ afd50d1c-3776-11eb-2696-abce7d038b88
begin
	x₀ = rand(3)
	n_steps = 2000
	s = [1.0, 4.0]
	orbit = run(solenoid, x₀, s, n_steps)
end

# ╔═╡ 3a605d5e-3776-11eb-34ff-67296e7637f3
let 
	plt = plot3d(
    1,
    xlim = (-3, 3),
    ylim = (-4, 3),
    zlim = (-0.7, 0.7),
    title = "Solenoid Attractor",
    marker = 2,
	leg = false,
	linealpha = 0
)
@gif for n = 1:1000
	push!(plt, orbit[n][1], orbit[n][2], orbit[n][3])
	end every 10
end

# ╔═╡ 5720aa0c-38ed-11eb-073e-2d5e0332329d
md"""
### Lyapunov exponents
"""

# ╔═╡ 7395d3e2-38ed-11eb-37e9-d734cb6fffae


# ╔═╡ Cell order:
# ╟─7abc4f04-3772-11eb-1c9b-712985ec2af7
# ╠═b49dd7c8-3772-11eb-3414-bb1f9e84b748
# ╠═f8e62384-3772-11eb-23fb-41462d88e988
# ╠═9d23a890-3773-11eb-1d5e-295158261eea
# ╠═afd50d1c-3776-11eb-2696-abce7d038b88
# ╠═3a605d5e-3776-11eb-34ff-67296e7637f3
# ╠═5720aa0c-38ed-11eb-073e-2d5e0332329d
# ╠═7395d3e2-38ed-11eb-37e9-d734cb6fffae
