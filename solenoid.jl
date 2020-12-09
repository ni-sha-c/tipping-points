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
	using Test
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
	
	r₁ = s₀ + (r - s₀)/s₁ + cos(θ)/2
	θ₁ = 2θ
	z₁ = z/s₁ + sin(θ)/2
	
	xn = similar(x)
	xn[1] = r₁*cos(θ₁)
	xn[2] = r₁*sin(θ₁)
	xn[3] = z₁
	
	return xn
end

# ╔═╡ 9d23a890-3773-11eb-1d5e-295158261eea
function run(method, x₀, s, n_steps)
	n = 1
	orbit = zeros(3, n_steps)
	x₁ = copy(x₀)
	while n <= n_steps
		x₁ = method(x₁, s)
		orbit[:,n] = x₁
		n += 1
	end
	return orbit
end

# ╔═╡ 86f96d6c-39b6-11eb-3e44-e3a3d68a6957
function spinup(method, x₀, s, n_steps)
	n = 1
	x₁ = similar(x₀)
	while n < n_steps
		x₁ = method(x₁, s)
		n += 1
	end
	return x₁
end

# ╔═╡ afd50d1c-3776-11eb-2696-abce7d038b88
begin
	x₀ = rand(3)
	n_steps = 1200
	s = [1.0, 4.0]
	x₀ = spinup(solenoid, x₀, s, 1000)
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
	linealpha = 0,
	xlabel = "x",
	ylabel = "y",
	zlabel = "z"
)
@gif for n = 1:n_steps
	push!(plt, orbit[1,n], orbit[2,n], orbit[3,n])
	end every 10
end

# ╔═╡ a4a93482-39b1-11eb-0753-c5aeb0eda8ae
let
	
	p1 = plot(orbit[1,:], orbit[2,:], m=:o, ms=2, linealpha=0, leg=false, xlabel="x", ylabel="y")
	p2 = plot(orbit[2,:], orbit[3,:], m=:o, ms=2, linealpha=0, leg=false, xlabel="y", ylabel="z")

	plot(p1, p2, layout = grid(1, 2, widths=[0.5, 0.5]))
end


# ╔═╡ 5720aa0c-38ed-11eb-073e-2d5e0332329d
md"""
### Lyapunov exponents
"""

# ╔═╡ c619c9e4-390a-11eb-3819-3d17ff4aecd2
function tangent(x, v, s=[1.,4.])
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
		
		
		v₁ = similar(v)
		v₁[1] = dx1_dr*dr_dx*v[1] + 
			   dx1_dθ*dθ_dx*v[1] + 
			   dx1_dr*dr_dy*v[2] + 
			   dx1_dθ*dθ_dy*v[2] 
			   
			   
		v₁[2] = dy1_dr*dr_dx*v[1] + 
			   dy1_dθ*dθ_dx*v[1] + 
			   dy1_dr*dr_dy*v[2] + 
			   dy1_dθ*dθ_dy*v[2] 
		
		v₁[3] = dz1_dθ*dθ_dx*v[1] +
			   dz1_dθ*dθ_dy*v[2] +
			   dz1_dz*v[3]
		return v₁
end

# ╔═╡ ca2dbc4e-3972-11eb-15e0-c595f6d2d087
function test_tangent(x, v, s=[1.,4.],ϵ=1.e-6)
	fx = solenoid(x, s)
	df_dx = solenoid(x .+ ϵ*[1., 0., 0.], s) - solenoid(x .- ϵ*[1., 0., 0.], s)
	df_dy = solenoid(x .+ ϵ*[0., 1., 0.], s) - solenoid(x .- ϵ*[0., 1., 0.], s)
	df_dz = solenoid(x .+ ϵ*[0., 0., 1.], s) - solenoid(x .- ϵ*[0., 0., 1.], s)
	v_ana = tangent(x, v)
	v_fd = Matrix([df_dx df_dy df_dz])*v/(2*ϵ)
	@show v_ana, v_fd
	@test v_ana ≈ v_fd atol=1.e-6
end

# ╔═╡ 7395d3e2-38ed-11eb-37e9-d734cb6fffae
function first_LE(x₀, s, n_steps)
	v = rand(3)
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

# ╔═╡ 83546b9c-39c2-11eb-0772-ff2a7d6ed0fe
function all_LEs(x₀, s, n_steps)
	v₁ = rand(3)
	v₂ = zeros(3)
	v₃ = zeros(3)
	Q = Matrix([v₁ v₂ v₃])
	R = similar(Q)
	x = copy(x₀)
	λ = zeros(3)
	n_spinup = 100
	α = 0.
	for n = 1:n_steps + n_spinup
		Q .= Matrix([v₁ v₂ v₃])
		A = qr(Q)
		Q .= Array(A.Q)
		R .= A.R
		@show R[1,1], norm(v₁)
		if n > n_spinup
			λ .+= log.(abs.(diag(R)))/n_steps
			α += log(norm(v₁))/n_steps
		end	
		
		v₁ .= Q[:,1]
		v₂ .= Q[:,2]
		v₃ .= Q[:,3]
		
		
		v₁ .= tangent(x, v₁)
		v₂ .= tangent(x, v₂)
		v₃ .= tangent(x, v₃)
		x .= solenoid(x, s)
	end
	@show α
	return λ
end

# ╔═╡ bbe580e0-3972-11eb-0fc7-4d95edb00ca0
test_tangent(rand(3), rand(3))

# ╔═╡ 204d842c-39ca-11eb-19fe-9fdffed8d7df
begin 
	n_le = 100
	x = spinup(solenoid, rand(3), s, n_le)
	all_LEs(x, s, 1000)
end

# ╔═╡ Cell order:
# ╟─7abc4f04-3772-11eb-1c9b-712985ec2af7
# ╠═b49dd7c8-3772-11eb-3414-bb1f9e84b748
# ╠═f8e62384-3772-11eb-23fb-41462d88e988
# ╠═9d23a890-3773-11eb-1d5e-295158261eea
# ╠═86f96d6c-39b6-11eb-3e44-e3a3d68a6957
# ╠═afd50d1c-3776-11eb-2696-abce7d038b88
# ╠═3a605d5e-3776-11eb-34ff-67296e7637f3
# ╠═a4a93482-39b1-11eb-0753-c5aeb0eda8ae
# ╠═5720aa0c-38ed-11eb-073e-2d5e0332329d
# ╠═c619c9e4-390a-11eb-3819-3d17ff4aecd2
# ╠═ca2dbc4e-3972-11eb-15e0-c595f6d2d087
# ╠═7395d3e2-38ed-11eb-37e9-d734cb6fffae
# ╠═83546b9c-39c2-11eb-0772-ff2a7d6ed0fe
# ╠═bbe580e0-3972-11eb-0fc7-4d95edb00ca0
# ╠═204d842c-39ca-11eb-19fe-9fdffed8d7df
