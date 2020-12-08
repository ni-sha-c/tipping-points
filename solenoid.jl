### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ b49dd7c8-3772-11eb-3414-bb1f9e84b748
begin
	using LaTeXStrings
	using Plots
	using PlutoUI
	using LinearAlgebra
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

# ╔═╡ c619c9e4-390a-11eb-3819-3d17ff4aecd2
function tangent(x, v)
		r_sq = x[1]*x[1] + x[2]*x[2]
		r = sqrt(r_sq)
		θ = (atan(x[2], x[1]) + 2π) % (2π) 
		z = x[3]
		
		ct, st = cos(θ), sin(θ)
		
		r₁ = s[1] + (r - s[1])/s[2] + ct/2
		θ₁ = 2θ
		z₁ = z/s[2] + st/2
		
		ct1, st1 = cos(θ₁), sin(θ₁)
		
		dx1_dr1 = ct1
		dx1_dθ1 = -r₁*st1
		dy1_dr1 = st1
		dy1_dθ1 = r₁*ct1
		
		dr1_dr = 1/s[2]
		dr1_dθ = -st/2
		dθ1_dθ = 2
		dz1_dz = 1/s[2]
		dz1_dθ = ct/2
		
		
		dr_dx = x[1]/r
		dr_dy = x[2]/r
		dθ_dx = -x[2]/r_sq
		dθ_dy = x[1]/r_sq
		
		dx1_dr = dx1_dr1*dr1_dr 
		dx1_dθ = dx1_dr1*dr1_dθ + dx1_dθ1*dθ1_dθ
		dy1_dr = dy1_dr1*dr1_dr 
		dy1_dθ = dy1_dr1*dr1_dθ + dy1_dθ1*dθ1_dθ
		
		
		
		v[1] = dx1_dr*dr_dx*v[1] + 
			   dx1_dθ*dθ_dx*v[1] + 
			   dx1_dr*dr_dy*v[2] + 
			   dx1_dθ*dθ_dy*v[2] 
			   
			   
		v[2] = dy1_dr*dr_dx*v[1] + 
			   dy1_dθ*dθ_dx*v[1] + 
			   dy1_dr*dr_dy*v[2] + 
			   dy1_dθ*dθ_dy*v[2] 
		
		v[3] = dz1_dθ*dθ_dx*v[1] +
			   dz1_dθ*dθ_dy*v[2]
			   dz1_dz*v[3]
		
		return v
end

# ╔═╡ 7395d3e2-38ed-11eb-37e9-d734cb6fffae
function LE(x₀, v₀, s, n_steps)
	v = copy(v₀)
	x = copy(x₀)
	λ = 0.
	for n = 1:n_steps
		v .= tangent(x, v)
		α = norm(v)
		λ += log(α)/n_steps
		v ./= α
		x .= solenoid(x, s)
	end
	return λ
end

# ╔═╡ 7d245e0e-390d-11eb-28dd-7b1504da3c0f
LE(rand(3), rand(3), [1.4, 1], 3500)

# ╔═╡ Cell order:
# ╟─7abc4f04-3772-11eb-1c9b-712985ec2af7
# ╠═b49dd7c8-3772-11eb-3414-bb1f9e84b748
# ╠═f8e62384-3772-11eb-23fb-41462d88e988
# ╠═9d23a890-3773-11eb-1d5e-295158261eea
# ╠═afd50d1c-3776-11eb-2696-abce7d038b88
# ╠═3a605d5e-3776-11eb-34ff-67296e7637f3
# ╠═5720aa0c-38ed-11eb-073e-2d5e0332329d
# ╠═c619c9e4-390a-11eb-3819-3d17ff4aecd2
# ╠═7395d3e2-38ed-11eb-37e9-d734cb6fffae
# ╠═7d245e0e-390d-11eb-28dd-7b1504da3c0f
