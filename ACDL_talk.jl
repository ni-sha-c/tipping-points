### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ f9250b62-3b5f-11eb-011a-ab7e41db50df
begin
	using LaTeXStrings
	using Plots
	using PlutoUI
	using LinearAlgebra
	using Test
end

# ╔═╡ f19a3c4c-3b48-11eb-1d00-fb92a48cfba2
md"""
### An introduction to studying climate phenomena on small systems
"""

# ╔═╡ a2ea8866-3b60-11eb-0001-ff9253168fdd
md"""
###### [Valerio Lucarini and Tamás Bódai. Global stability properties of the climate: Melancholia states, invariant measures, and phase transitions Nonlinearity 33 R59 (2020) ](http://dx.doi.org/10.1088/1361-6544/ab86cc)

###### [An accessible explainer on climate tipping points](https://www.carbonbrief.org/explainer-nine-tipping-points-that-could-be-triggered-by-climate-change)

###### [Kiers, C., Jones, C.K.R.T. On Conditions for Rate-induced Tipping in Multi-dimensional Dynamical Systems. J Dyn Diff Equat 32, 483–503 (2020)](https://doi.org/10.1007/s10884-019-09730-9)

###### [Ashwin, P., Wieczorek, S., Vitolo, R., Cox, P.: Tipping points in open systems: bifurcation, noise-induced and rate-dependent examples in the climate system. Philos. Trans. R. Soc. Lond. A: Math. Phys. Eng.Sci.370, 1166–1184 (2012)](https://doi.org/10.1098/rsta.2011.0306)

###### [Lenton, T. Early warning of climate tipping points. Nature Climate Change Vol 1 (2011)](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.666.244&rep=rep1&type=pdf)

"""



# ╔═╡ f0332982-3b62-11eb-34d9-9fbf204d3c99
md"""

##### "A climate ‘tipping point’ occurs when a small change in forcing triggers a strongly nonlinear response in the internal dynamics of part of the climate system, qualitatively changing its future state." -- Lenton, T. 2011 

"""

# ╔═╡ 4f41e004-3b62-11eb-29f3-497135307b6a
md"""
In [Ashwin et al.'s](https://doi.org/10.1098/rsta.2011.0306) classification, bifurcation-induced tipping points, and more recently noise-induced and rate-induced tipping points have been studied for early warning signs.
These studies typically involve the computation of critical stability exponents of trajectories. 

"""

# ╔═╡ aec084c4-3b67-11eb-2cbe-c387e6480374
md"""
#### Parametric dependence of the Solenoid attractor
"""

# ╔═╡ 2d5230d2-3b5f-11eb-1d2b-61caca80e341
md"""
##### The Smale-Williams Solenoid map is a three-dimensional chaotic map that produces a *uniformly hyperbolic* attractor. Uniform hyperbolicity is the simplest mathematical setting in which chaotic attractors occur. 

Uniform hyperbolicity refers to the presence of uniformly (exponentially) expanding and contracting tangent subspaces at every phase point.
"""

# ╔═╡ e855fa16-3b5e-11eb-0e39-f1afae7c0a09
function random_solenoid(x, s, p=0.5)
	
	s₀, s₁ = s[1], s[2]
	if rand() < p 
		s₀ = 2s₀
	end
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

# ╔═╡ f691f898-3b5e-11eb-3b9c-bf371cbe226f
function solenoid(x, s)
	return random_solenoid(x, s, 0)
end

# ╔═╡ 941de504-3b5f-11eb-1e98-c52dda487540
md"""
### Function Library

"""

# ╔═╡ 0c378bdc-3b60-11eb-371a-89bb7f6d480f
md"""
#### Dependencies
"""

# ╔═╡ a65139e2-3b5f-11eb-1cd2-f34cf61bab4c
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

# ╔═╡ b00d7fae-3b5f-11eb-3ff6-e1384b3ceb8c
function run(method, x₀, s, n_steps, p)
	n = 1
	orbit = zeros(3, n_steps)
	x₁ = copy(x₀)
	while n <= n_steps
		x₁ = method(x₁, s, p)
		orbit[:,n] = x₁
		n += 1
	end
	return orbit
end

# ╔═╡ b62d21f0-3b5f-11eb-1633-7df578a7f45b
function spinup(method, x₀, s, n_steps)
	n = 1
	x₁ = copy(x₀)
	while n < n_steps
		x₁ = method(x₁, s)
		n += 1
	end
	return x₁
end

# ╔═╡ 2a46f580-3b5f-11eb-02b5-5b0affe37c11
begin
	x₀ = rand(3)
	n_steps = 5000
	s = [1.0, 4.0]
	x₀ = spinup(random_solenoid, x₀, s, 1000)
	orbit_r = run(random_solenoid, x₀, s, n_steps)
	orbit = run(solenoid, x₀, s, n_steps)
end

# ╔═╡ d0c62c46-3b5f-11eb-3df6-4165a5d3db41
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

# ╔═╡ b0d93120-3b65-11eb-244d-db50471180a6
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
	push!(plt, orbit_r[1,n], orbit_r[2,n], orbit_r[3,n])
	end every 10
end

# ╔═╡ e2e3875e-3b67-11eb-2227-cb921c13814c
let
	
	p1 = plot(orbit[1,:], orbit[2,:], m=:o, ms=2, linealpha=0, leg=false, xlabel="x", ylabel="y")
	p2 = plot(orbit[2,:], orbit[3,:], m=:o, ms=2, linealpha=0, leg=false, xlabel="y", ylabel="z")

	plot(p1, p2, layout = grid(1, 2, widths=[0.5, 0.5]))
end


# ╔═╡ fadc9fee-3b67-11eb-0ff2-6ff1b81fd678
let
	
	p1 = plot(orbit_r[1,:], orbit_r[2,:], m=:o, ms=2, linealpha=0, leg=false, xlabel="x", ylabel="y")
	p2 = plot(orbit_r[2,:], orbit_r[3,:], m=:o, ms=2, linealpha=0, leg=false, xlabel="y", ylabel="z")

	plot(p1, p2, layout = grid(1, 2, widths=[0.5, 0.5]))
end


# ╔═╡ Cell order:
# ╟─f19a3c4c-3b48-11eb-1d00-fb92a48cfba2
# ╟─a2ea8866-3b60-11eb-0001-ff9253168fdd
# ╟─f0332982-3b62-11eb-34d9-9fbf204d3c99
# ╟─4f41e004-3b62-11eb-29f3-497135307b6a
# ╟─aec084c4-3b67-11eb-2cbe-c387e6480374
# ╟─2d5230d2-3b5f-11eb-1d2b-61caca80e341
# ╟─e855fa16-3b5e-11eb-0e39-f1afae7c0a09
# ╟─f691f898-3b5e-11eb-3b9c-bf371cbe226f
# ╠═2a46f580-3b5f-11eb-02b5-5b0affe37c11
# ╟─d0c62c46-3b5f-11eb-3df6-4165a5d3db41
# ╟─b0d93120-3b65-11eb-244d-db50471180a6
# ╟─e2e3875e-3b67-11eb-2227-cb921c13814c
# ╟─fadc9fee-3b67-11eb-0ff2-6ff1b81fd678
# ╟─941de504-3b5f-11eb-1e98-c52dda487540
# ╟─0c378bdc-3b60-11eb-371a-89bb7f6d480f
# ╟─f9250b62-3b5f-11eb-011a-ab7e41db50df
# ╟─a65139e2-3b5f-11eb-1cd2-f34cf61bab4c
# ╟─b00d7fae-3b5f-11eb-3ff6-e1384b3ceb8c
# ╟─b62d21f0-3b5f-11eb-1633-7df578a7f45b
