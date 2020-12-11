### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

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

# ╔═╡ f3a666d6-3b69-11eb-1b15-7bada2853192
md"""
We denote our chaotic maps $\varphi_s$, where $s$ is the set of parameters.
The maps have a compact hyperbolic set and hence an SRB measure over a range of $s$ values.

$\varphi_s^n$ refers to the $n$-time composition of $\varphi_s$. That is, $\varphi^n_s = \varphi_s \circ \varphi^{n-1}_s$.

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

# ╔═╡ ca68ef88-3b68-11eb-1dc4-11574ee92984
md"""
### Lyapunov Exponents
"""

# ╔═╡ c8bff564-3b68-11eb-134a-598082ea5f4f
md"""
Lyapunov epxponents measure the asymptotic sensitivity to infinitesimal perturbations.

The tangent space splits as a direct sum $\oplus E_i(x)$, at every $x$, such that  tangent vectors in $E_i$ are covariant and have the same asymptotic exponential growth rate, which are called Lyapunov exponents.
$\lambda_i  = \lim_{N\to \infty} \dfrac{1}{N} \log \| d\varphi^n_s (x) E_i(x) \|$
"""

# ╔═╡ aa67bff0-3b6a-11eb-31ad-9decdab06358
md""" 
### First parameter variation
"""

# ╔═╡ 78911fc0-3b6b-11eb-1d24-0dfa6c85c165
md""" 
### Second parameter variation
"""

# ╔═╡ 578c1748-3b6c-11eb-0e4e-c5f2b86d0558
md"""
LEs in this case are not indicative of marked changes in the attractor.
"""

# ╔═╡ 8225f49c-3b6c-11eb-334f-056b4a1e4d4d
md""" 
### Frobenius-Perron Approach: a measure of global parametric dependence? 
"""

# ╔═╡ b30bdba8-3b6c-11eb-3551-39397136cdb7
md"""
The Frobenius-Perron operator $P:L^1\to L^1$ evolves probability densities under $\varphi$. For any set $A$, and a bounded function $g$,
 $\int_A g\; P f \; dx = \int_{\varphi^{-1} A} f\; (g\circ\varphi)\; dx.$
This operator lets us analyze ensemble behavior.
"""

# ╔═╡ bdcc30ce-3b6c-11eb-064a-0791a37f9bbd
md"""
In hyperbolic systems, $P$ may not be compact. But, we nevertheless try to approximate this operator on a finite-dimensional indicator function basis. Numerically this can be done by dividing up the domain into smaller cells and constructing the Markov transition matrix 

"""

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


# ╔═╡ 5685d37e-3b68-11eb-2abd-b5e245911928
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

# ╔═╡ 7baf192e-3b68-11eb-0b71-a9987e96a0df
function all_random_LEs(x₀, s, n_steps,p=0.5)
	v₁ = rand(3)
	v₂ = zeros(3)
	v₃ = zeros(3)
	Q = Matrix([v₁ v₂ v₃])
	R = similar(Q)
	x = copy(x₀)
	λ = zeros(3)
	n_spinup = 100
	α = 0.
	s₀, s₁ = s
	for n = 1:n_steps + n_spinup
		Q .= Matrix([v₁ v₂ v₃])
		A = qr(Q)
		Q .= Array(A.Q)
		R .= A.R
		if n > n_spinup
			λ .+= log.(abs.(diag(R)))/n_steps
			α += log(norm(v₁))/n_steps
		end	
		
		v₁ .= Q[:,1]
		v₂ .= Q[:,2]
		v₃ .= Q[:,3]
		
		
		v₁ .= tangent(x, v₁, [s₀, s₁])
		v₂ .= tangent(x, v₂, [s₀, s₁])
		v₃ .= tangent(x, v₃, [s₀, s₁])
		s₀ = s[1]
		if rand() < p
			s₀ = 2s₀
		end
		x .= solenoid(x, [s₀, s₁])
	end
	@show α
	return λ
end

# ╔═╡ 871c0df8-3b68-11eb-3431-578143a5007b
function all_LEs(x₀, s, n_steps)
	return all_random_LEs(x₀, s, n_steps, 0)
end

# ╔═╡ 4a247006-3b6b-11eb-3c7d-d16cea1ba192
function compute_LEs(s) 
	n_le = 100
	x = spinup(solenoid, rand(3), s, n_le)
	return all_LEs(x, s, 1000)
end

# ╔═╡ debf5cd6-3b6a-11eb-1496-7b9e0ab7e406
begin
	n_p = 10
	les_1 = zeros(3, n_p)
	s1 = LinRange(1.0, 6.5, n_p)
	for i = 1:n_p
		les_1[:,i] = compute_LEs([s1[i], 4.0])
	end
	les_1 = les_1'
	p1 = plot(s1, les_1[:,1], m=:o, title="1st LE", leg=false)
	p2 = plot(s1, les_1[:,2], m=:o, xlim=(1.0, 2.6), ylim=(-2, -1), title="2nd LE", leg=false)
	p3 = plot(s1, les_1[:,3], m=:o, xlim=(1.0, 2.6), ylim=(-2, -1), title="3rd LE", leg=false)
	
	plot(p1, p2, p3, layout = grid(1, 3, widths=[0.3, 0.3, 0.3]))
end

# ╔═╡ b9e78fc2-3b6b-11eb-2023-b138c8257cf2
@bind n Slider(1:n_p, show_value=true)

# ╔═╡ fd84434c-3b6b-11eb-2b6b-fbe51b601ea6
@bind m Slider(1:n_p, show_value=true)

# ╔═╡ 264bb3a0-3b6c-11eb-1d2d-476662b86ced
begin
	p10 = plot(s1[1:m], les_1[1:m,1], m=:o, xlim=(minimum(s1), maximum(s1)), ylim=(minimum(les_1[:,1]), maximum(les_1[:,1])),title="1st LE", xlabel=L"s_1", ylabel=L"\lambda_1", leg=false)
end

# ╔═╡ 3e1ccf50-3b6c-11eb-04ec-b51fba6a64c9
begin
	x₀ .= spinup(solenoid, x₀, [s1[m], 4.0], 5000)
	orbit_s1 = zeros(n_steps, 3, n_p)
	orbit_s1[:,:,m] = run(solenoid, x₀, [s1[m], 4.0], n_steps)'
	p11 = plot(orbit_s1[:,1,m], orbit_s1[:,2,m], m=:o, ms=2, linealpha=0, xlim=(-7.5,7.5), ylim=(-7.5,7.5),leg=false, xlabel="x", ylabel="y")
	p12 = plot(orbit_s1[:,2,m], orbit_s1[:,3,m], m=:o, ms=2, linealpha=0,
		xlim=(-7.5,7.5), ylim=(-1.,1.),leg=false, xlabel="y", ylabel="z")

	plot(p11, p12, layout = grid(1, 2, widths=[0.5, 0.5]))
end

# ╔═╡ 4ae256aa-3b6d-11eb-19fb-251716c33c79
@bind j Slider(1:n_p, show_value=true)

# ╔═╡ 9895cece-3b6b-11eb-21ee-af0764651cc9
begin
	s2 = LinRange(2.0, 10., n_p)
	les_2 = zeros(n_p, 3)
	for i = 1:n_p
		les_2[i,:] = compute_LEs([1.0, s2[i]])
	end
	
	p4 = plot(s2, les_2[:,1], m=:o,title="1st LE", leg=false)
	p5 = plot(s2, les_2[:,2], m=:o, title="2nd LE", leg=false)
	p6 = plot(s2, les_2[:,3], m=:o, title="3rd LE", leg=false)
	
	plot(p4, p5, p6, layout = grid(1, 3, widths=[0.3, 0.3, 0.3]))
end

# ╔═╡ d3c330ae-3b6b-11eb-385f-1d5068855986
begin
	p7 = plot(s2[1:n], les_2[1:n,1], m=:o, ylim=(minimum(les_2[:,1]), maximum(les_2[:,1])), xlim=(minimum(s2), maximum(s2)), title="1st LE", leg=false)
end

# ╔═╡ e87c0494-3b6b-11eb-2234-abcc80a22b2a
begin
	x₀ .= spinup(solenoid, x₀, [1.0, s2[n]], 500)
	orbit_s2 = zeros(n_steps, 3, n_p)
	orbit_s2[:,:,n] = run(solenoid, x₀, [1.0, s2[n]], n_steps)'
	p8 = plot(orbit_s2[:,1,n], orbit_s2[:,2,n], m=:o, ms=2, linealpha=0, leg=false, xlabel="x", ylabel="y")
	p9 = plot(orbit_s2[:,2,n], orbit_s2[:,3,n], m=:o, ms=2, linealpha=0, leg=false, xlabel="y", ylabel="z")

	plot(p8, p9, layout = grid(1, 2, widths=[0.5, 0.5]))
end

# ╔═╡ 56233a2c-3b6b-11eb-2deb-cdad6683a639
function compute_random_LEs(s,p=0.5) 
	n_le = 100
	x = spinup(solenoid, rand(3), s, n_le)
	return all_random_LEs(x, s, 3000, p)
end

# ╔═╡ 6e265c80-3b6b-11eb-1a9d-05a92e38c045
begin
	les_r1 = zeros(3, n_p)
	for i = 1:n_p
		les_r1[:,i] = compute_random_LEs([s1[i], 4.0],0.9)
	end
	les_r1 = les_r1'
	p14 = plot(s1, les_r1[:,1], m=:o, title="1st LE", leg=false)
	p15 = plot(s1, les_r1[:,2], m=:o, title="2nd LE", leg=false)
	p16 = plot(s1, les_r1[:,3], m=:o, title="3rd LE", leg=false)
	
	plot(p14, p15, p16, layout = grid(1, 3, widths=[0.3, 0.3, 0.3]))
end

# ╔═╡ b04f349c-3b6b-11eb-3539-c5dd7bbd7a72
begin
	
	les_r2 = zeros(n_p, 3)
	for i = 1:n_p
		les_r2[i,:] = compute_random_LEs([1.0, s2[i]])
	end
	
	p17 = plot(s2, les_r2[:,1], m=:o,title="1st LE", leg=false)
	p18 = plot(s2, les_r2[:,2], m=:o, title="2nd LE", leg=false)
	p19 = plot(s2, les_r2[:,3], m=:o, title="3rd LE", leg=false)
	
	plot(p17, p18, p19, layout = grid(1, 3, widths=[0.3, 0.3, 0.3]))
end

# ╔═╡ c7493b10-3b6c-11eb-21ce-274b1436b801
function construct_transition_matrix(orbit, n_nodes)
	P = zeros(n_nodes^3, n_nodes^3)
	n = size(orbit)[2]
	eps = 0.01
	x_min, x_max = minimum(orbit[1,:]) - eps, maximum(orbit[1,:]) + eps
	y_min, y_max = minimum(orbit[2,:]) - eps, maximum(orbit[2,:]) + eps
	z_min, z_max = minimum(orbit[3,:]) - eps, maximum(orbit[3,:]) + eps
	dx = (x_max - x_min)/n_nodes
	dy = (y_max - y_min)/n_nodes
	dz = (z_max - z_min)/n_nodes
	
	
	x, y, z = orbit[:,1]
	pre_ind_x = ceil(Int64,(x - x_min)/dx)
	pre_ind_y = ceil(Int64,(y - y_min)/dy)
	pre_ind_z = ceil(Int64,(z - z_min)/dz)
	pre_ind = pre_ind_x + (pre_ind_y-1)*n_nodes + (pre_ind_z-1)*n_nodes*n_nodes
	
	for i = 2:n
		
		x, y, z = orbit[:,i]
		
		ind_x = ceil(Int64,(x - x_min)/dx)
		ind_y = ceil(Int64,(y - y_min)/dy)
		ind_z = ceil(Int64,(z - z_min)/dz)
		
		ind = ind_x + (ind_y-1)*n_nodes + (ind_z-1)*n_nodes*n_nodes
		
		
		P[ind, pre_ind] += 1
		pre_ind = ind
	end
	for i = 1:n_nodes^3
		a = sum(P[:,i])
		if a > 0
			P[:, i] ./= sum(P[:,i]) 
		end
	end
	return P'
end

# ╔═╡ e28bf2e6-3b6c-11eb-1ce2-51345710e0ef
begin
	n_nodes = 5
	P1_r = zeros(n_nodes^3, n_nodes^3, n_p)
	for (i, s1i) in enumerate(s1)
		@show s1i
		x = spinup(solenoid, rand(3), [s1i, 4.0], 1000)
		test_orbit = run(random_solenoid, x, [s1i, 4.0], 2000000)
		P1_r[:,:,i]  = construct_transition_matrix(test_orbit, n_nodes)
	end
end

# ╔═╡ f1a7adc4-3b6c-11eb-1588-eb2a5c6a400a
begin
	rpr1_r = 1im*zeros(n_nodes^3, n_p)
	for i =1:n_p
		rpr1_r[:,i] = eigvals(P1_r[:,:,i])
	end
end

# ╔═╡ 35df2fa6-3b6d-11eb-210a-c5146dbec4b4
begin
	p20 = plot(real(rpr1_r[:,j]), imag(rpr1_r[:,j]), xlim=(-1,1), ylim=(-1,1), linealpha=0, ms=2, m=:o, leg=false, aspect_ratio=1)
	t_arr = 0.:pi/20:2π
	plot!(p20, cos.(t_arr), sin.(t_arr))
end

# ╔═╡ Cell order:
# ╟─f19a3c4c-3b48-11eb-1d00-fb92a48cfba2
# ╟─a2ea8866-3b60-11eb-0001-ff9253168fdd
# ╟─f0332982-3b62-11eb-34d9-9fbf204d3c99
# ╟─4f41e004-3b62-11eb-29f3-497135307b6a
# ╟─aec084c4-3b67-11eb-2cbe-c387e6480374
# ╟─2d5230d2-3b5f-11eb-1d2b-61caca80e341
# ╟─f3a666d6-3b69-11eb-1b15-7bada2853192
# ╟─e855fa16-3b5e-11eb-0e39-f1afae7c0a09
# ╟─f691f898-3b5e-11eb-3b9c-bf371cbe226f
# ╠═2a46f580-3b5f-11eb-02b5-5b0affe37c11
# ╟─d0c62c46-3b5f-11eb-3df6-4165a5d3db41
# ╟─b0d93120-3b65-11eb-244d-db50471180a6
# ╟─e2e3875e-3b67-11eb-2227-cb921c13814c
# ╟─fadc9fee-3b67-11eb-0ff2-6ff1b81fd678
# ╠═ca68ef88-3b68-11eb-1dc4-11574ee92984
# ╠═c8bff564-3b68-11eb-134a-598082ea5f4f
# ╠═aa67bff0-3b6a-11eb-31ad-9decdab06358
# ╟─debf5cd6-3b6a-11eb-1496-7b9e0ab7e406
# ╟─6e265c80-3b6b-11eb-1a9d-05a92e38c045
# ╟─78911fc0-3b6b-11eb-1d24-0dfa6c85c165
# ╟─9895cece-3b6b-11eb-21ee-af0764651cc9
# ╟─b04f349c-3b6b-11eb-3539-c5dd7bbd7a72
# ╟─b9e78fc2-3b6b-11eb-2023-b138c8257cf2
# ╠═d3c330ae-3b6b-11eb-385f-1d5068855986
# ╠═e87c0494-3b6b-11eb-2234-abcc80a22b2a
# ╠═fd84434c-3b6b-11eb-2b6b-fbe51b601ea6
# ╠═264bb3a0-3b6c-11eb-1d2d-476662b86ced
# ╟─3e1ccf50-3b6c-11eb-04ec-b51fba6a64c9
# ╟─578c1748-3b6c-11eb-0e4e-c5f2b86d0558
# ╠═8225f49c-3b6c-11eb-334f-056b4a1e4d4d
# ╠═b30bdba8-3b6c-11eb-3551-39397136cdb7
# ╠═bdcc30ce-3b6c-11eb-064a-0791a37f9bbd
# ╠═e28bf2e6-3b6c-11eb-1ce2-51345710e0ef
# ╠═f1a7adc4-3b6c-11eb-1588-eb2a5c6a400a
# ╠═4ae256aa-3b6d-11eb-19fb-251716c33c79
# ╠═35df2fa6-3b6d-11eb-210a-c5146dbec4b4
# ╟─941de504-3b5f-11eb-1e98-c52dda487540
# ╟─0c378bdc-3b60-11eb-371a-89bb7f6d480f
# ╟─f9250b62-3b5f-11eb-011a-ab7e41db50df
# ╟─a65139e2-3b5f-11eb-1cd2-f34cf61bab4c
# ╟─b00d7fae-3b5f-11eb-3ff6-e1384b3ceb8c
# ╟─b62d21f0-3b5f-11eb-1633-7df578a7f45b
# ╟─5685d37e-3b68-11eb-2abd-b5e245911928
# ╟─7baf192e-3b68-11eb-0b71-a9987e96a0df
# ╟─871c0df8-3b68-11eb-3431-578143a5007b
# ╟─4a247006-3b6b-11eb-3c7d-d16cea1ba192
# ╟─56233a2c-3b6b-11eb-2deb-cdad6683a639
# ╟─c7493b10-3b6c-11eb-21ce-274b1436b801
