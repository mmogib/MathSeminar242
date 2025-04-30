### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 59a80010-21a8-11f0-2aaa-5528a16a7081
begin
	
	using Pkg: Pkg
	Pkg.activate("../")
	using Plots, PlotThemes
	using OAFCVX
	using CommonMark
	using PlutoUI
	using LaTeXStrings, Latexify, HypertextLiteral, Colors
	using LinearAlgebra, Dates
	using JuMP, Ipopt
	md"``\color{white}{\text{Packages}}``"
end

# ╔═╡ 1c7f4f9c-0954-408f-927f-66d5929fc006
TableOfContents(title="📚 MATH SEMINAR: TERM 242", indent=true, depth=4)

# ╔═╡ 85b8d34e-d89d-49f0-b44e-fc0acee01d73
md"# Approximating Convex Vector Optimization Solutions"

# ╔═╡ 8c4c7c31-e590-4ff5-bc7b-f14d15293980
cm"""
This seminar provides an accessible introduction to scalarization techniques for approximating solutions of convex vector optimization problems. After motivating the general convex vector‐optimization framework and its “upper image,” we review key definitions—C-convexity, minimal and weakly minimal elements, polyhedral cones, and ε-solutions—through intuitive examples. Emphasis is placed on clear geometric interpretation and practical implementation. We conclude by discussing promising research directions: handling unbounded vector problems and establishing precise rates of convergence for these algorithms.
"""

# ╔═╡ dc6dbd80-444f-4697-ad54-476b16cd5fe6
md"## The Problem"

# ╔═╡ c84c7d54-89b8-481e-9b2b-26ede92942f6
md"## Example"

# ╔═╡ 764259e7-1d76-495b-8553-89f72caf3e46
cm"## Definitions and Background"

# ╔═╡ 1c060d79-1857-495b-994a-90273540afde
md"## Solutions Methods"

# ╔═╡ fcb71315-2400-4c64-89ed-a68de52c2cbd
md"### Weighted Sum"

# ╔═╡ d17c513c-768d-4dba-98ae-e7312e84535b
cm"""
For a weight parameter ``w \in \mathbb{R}^m``, the __weighted sum scalarization model__ is given by
```math
\begin{array}{ll}
\min & w^{\top} f(x) \\
\text { subject to } & x \in \mathcal{\Omega}.
\end{array}\tag{WSw}
```
"""

# ╔═╡ d13e1a2d-bb95-4a45-b920-41abe2ad1a3e
md"### Pascoletti–Serafini"

# ╔═╡ 204efc5c-9834-4e28-ac2b-26bd1a659a1c
cm"""
For ``v, d \in \mathbb{R}^m``, __Pascoletti–Serafini scalarization__ is given by
```math
\begin{array}{ll}
\min & z \\
\text { subject to } \\
& f(x) \leq_C v + zd\\
& x \in \mathcal{\Omega}\\
& z \in \mathbb{R}\\
\end{array}\tag{PS(v,d)}
```

In what follows, we need the Lagrange dual of (PS(v,d)) which is give by
```math
\begin{array}{ll}
\max & \displaystyle\inf_{x\in \Omega}w^T(f(x)-v) \\
\text { subject to } \\
& w^T d = 1\\
& w \in C^+\\
\end{array}\tag{DPS(v,d)}
```
"""

# ╔═╡ f9a1dfe1-bc71-42c4-a6e8-3f422eba7507
# let
# 	#── 1) draw P = [0,2.5]^2 ──#
# 	rect = Shape([0.0, 2.5, 2.5, 0.0], [0.0, 0.0, 2.5, 2.5])
# 	p = plot(rect;
# 	    seriestype   = :shape,
# 	    fillcolor    = :lightgray,
# 	    fillalpha    = 0.6,
# 	    linecolor    = :transparent,
# 	    xlims        = (-0.2, 2.6),
# 	    ylims        = (-0.2, 2.6),
# 	    aspect_ratio = :equal,
# 	    xlabel       = L"f_1",
# 	    ylabel       = L"f_2",
# 	    legend       = false,
# 	    grid         = false,
# 		frame_style=:origin
# 	)
	
# 	#── 2) overlay the disk f(X) ──#
# 	xc, yc, r = 1.0, 1.0, 1.0
# 	θ = range(0, 2π, length=300)
# 	disk = Shape(xc .+ r*cos.(θ), yc .+ r*sin.(θ))
# 	plot!(p, disk;
# 	    seriestype   = :shape,
# 	    fillcolor    = :gray,
# 	    fillalpha    = 0.8,
# 	    linecolor    = :darkgray,
# 	    linewidth    = 1.5,
# 	)
	
# 	#── 3) pick the tangency point p on the circle at 225° ──#
# 	θ_t = 5π/4
# 	p_x = xc + r*cos(θ_t)
# 	p_y = yc + r*sin(θ_t)
# 	scatter!(p, [p_x], [p_y];
# 	    marker      = :circle,
# 	    color       = :blue,
# 	    markersize  = 2,
# 	)
	
# 	#── 4) the inward normal d (unit) ──#
# 	d = -[cos(θ_t), sin(θ_t)]  # points toward the center
# 	quiver!(p, [0.0], [0.0],
# 	    quiver      = ([d[1]], [d[2]]),
# 	    arrow       = true,
# 	    linewidth   = 1.2,
# 	    linecolor   = :blue,
# 	    linestyle   = :dot,
# 	)
	
# 	#── 5) hyperplane H: { x | d⋅x = d⋅p } ──#
# 	# compute intercept
# 	b = dot(d, [p_x, p_y])
# 	xs = [-0.2, 2.6]
# 	ys = (b .- d[1]*xs) ./ d[2]
# 	plot!(p, xs, ys;
# 	    linecolor  = :red,
# 	    linestyle  = :dash,
# 	    linewidth  = 2,
# 	)
	
# 	#── 6) support‐function points on axes ──#
# 	# intercept at x‐axis: (b/d₁, 0), at y‐axis: (0, b/d₂)
# 	pt_x, pt_y = b/d[1], b/d[2]
# 	scatter!(p, [pt_x, 0.0], [0.0, pt_y];
# 	    color      = :red,
# 	    markersize = 6,
# 	)
	
# 	#── 7) mark the origin v = 0 ──#
# 	scatter!(p, [0.0], [0.0];
# 	    color      = :black,
# 	    markersize = 4,
# 	)
	
# 	#── 8) add annotations ──#
# 	annotate!(p, p_x+0.1,  p_y+0.1,  text("p",               :blue, 12))
# 	annotate!(p, p_x+0.4,  p_y+0.4,  text(L"d",             :blue, 14))
# 	annotate!(p, -0.1,     -0.1,     text(L"v=0",            :black,12))
# 	annotate!(p, 0.0,      pt_y+0.2, text(L"\mathcal{H}",   :red,  14))
# 	annotate!(p, 1.8,      2.2,      text(L"\mathcal{P}",   :black,16))	
# end

# ╔═╡ e3f5835c-18fc-4f00-92f5-e55598e68081
md"## The algorithm"

# ╔═╡ c6db5f72-a8fa-481c-85e9-7c587e0589b4
cm"""
### Initialization
<div style="font-size:1.5em;line-height:2.5em;">

- Let ``w_i, 1\leq i\leq l`` be __extreme directions__ of ``C^+`` and
```math
x_i = \operatorname{solution}\; \text{WSw}_i, \quad \text{for }i \in \{1,2,\cdots, l\}.
```
- Let ``\displaystyle \overline{\Omega}_0=\{x_1,\cdots,x_l\}``.
- Compute the initial outer approximation ``P_0`` of ``\mathcal{P}`` as follows
```math
P_0 = \bigcap_{i=1}^l \mathcal{H}_i = \bigcap_{i=1}^l \left\{y\in \mathbb{R}^m\;|\; w_i^Ty\geq w_i^Tf(x_i)\right\} 
```
- Compute the set of vertices ``V_0`` of ``P_0``.
- Let ``V_{\text{used}}=\emptyset`` be the set of used (visited) vertices.

### Main Loop (k):
Let ``V_k`` be the set of vertices of the current outer approximation ``P_k``.

__While ``V_k \setminus V_{\text{used}}\not=\emptyset`` do__

- __Choose an unused vertex__ ``v \in V_k \setminus V_{\text{used}}``.
- __Choose a direction__ ``d_v \in \operatorname{int} C`` 
- __Solve ``PS(v,d_v)`` and ``DPS(v,d_v)``__ and
```math
\begin{array}{rcl}
x_v, z_v &=& \operatorname{solution pair}\; \text{PS(v,d}_v\text{)}.\\
w_v &=& \operatorname{solution}\; \text{D-PS(c,d}_v\text{)}.\\
\end{array}
```
- Update: ``V_{\text{used}}=V_{\text{used}}\bigcup \{v\}, \quad \overline{\Omega}_k=\overline{\Omega}_k \bigcup \{x_v\}``. 
- If ``z_v> \epsilon`` then
  - Compute ``\mathcal{H} = \left\{y\in \mathbb{R}^m\;|\; w_v^Ty\geq w_v^Tv+z_v\right\}``
  - Update ``P_{k+1} = P_k \bigcap \mathcal{H},\quad \overline{\Omega}_{k+1}=\overline{\Omega}_{k}``
  - Compute the set of vertices ``V_{k+1}`` of ``P_{k+1}``.
  - ``k = k + 1``
</div>
"""

# ╔═╡ 574605cb-46c1-4470-b4aa-b4a5e6fe6f8a
md"""
## Solving the Example 
"""

# ╔═╡ 456f82ae-67de-4050-8d47-e1c106bdd5d1
cm"""
__Initialization__
- Let ``w_1=(1,0), w_2=(0,1)`` be __extreme directions__ of ``C^+=\mathbb{R}^2_+`` and
```math
x_1 = (1,2), \quad x_2 = (2,1)
```
- Let ``\displaystyle \overline{\Omega}_0=\{x_1,x_2\}``.
- Compute the initial outer approximation ``P_0`` of ``\mathcal{P}`` as follows
```math
\begin{array}{lcl}
P_0 &=& \bigcap_{i=1}^2 \mathcal{H}_i = \bigcap_{i=1}^2 \left\{y\in \mathbb{R}^2\;|\; w_i^Ty\geq w_i^Tf(x_i)\right\} \\
&=&\{(x,y)\in \mathbb{R}^2\;|\;x\geq 1\}\bigcap\{(x,y)\in \mathbb{R}^2\;|\;y\geq 1\}\\
\end{array}
```
- Compute the set of vertices ``V_0`` of ``P_0``.
```math
V_0=\{(1,1)\}
```
- Let ``V_{\text{used}}=\emptyset`` be the set of used (visited) vertices.

"""

# ╔═╡ 22e078de-0f08-4dd9-bb15-333f5cd7f22b
cm"__Run Algorithm__ step by step"

# ╔═╡ b5d6e08e-318b-4f30-96c3-fe096f6da8dd
begin
	steps_ = @bind steps Radio(["1"=>"Initial", "2"=>"Iteration 1", "3"=>"Iteration 5"], default="1")
	# cm"""
	# __Run Algorithm__ step by step

	# $(steps_)
	# """
end

# ╔═╡ 453c1441-2b72-4650-a317-56aec0def692
md"## Choosing A Direction"

# ╔═╡ f82ea70e-0c04-4275-afbf-4cff107940bf
cm"""
1. __Fixed Direction__.
1. __Adjacent vertices-based approach__
2. __Ideal point-based approach__
"""

# ╔═╡ 44580f32-df76-484b-9d2f-d09723747e66
md"## Choosing A Vertex"

# ╔═╡ a9a7ac5b-6f09-4c90-a867-aee556f27c26
cm"""
1. __Vertex selection with clusters__
2. __Vertex selection with adjacency information__
3. __Vertex selection using local upper bounds__
"""

# ╔═╡ ef9f9816-74df-4d2d-ac27-12ccd4c1c755
md"## Convergence analysis"

# ╔═╡ 513d75b6-63ea-40f4-89b1-f1c79ac83c88
cm"""
If (PS(v,d)) is replaced by the __norm minimization scalarization__
```math
\begin{array}{ll}
\min & \|z\| \\
\text { subject to } \\
& f(x) \leq_C v + z\\
& x \in \mathcal{\Omega}\\
& z \in \mathbb{R}^m\\
\end{array}\tag{NM(v)}
```

> Ararat et. al. in __Convergence Analysis of a Norm Minimization-Based Convex Vector Optimization Algorithm. (n.d.). https://doi.org/10.1137/23M1574580 (2024)__
proved that for an arbitrary norm used in the scalarization models, 
> - the approximation error after ``k`` iterations decreases by the order of ``\mathcal{O}\left(k^{1 /(1-m)}\right)``, where ``m`` is the dimension of the objective space. 
> - in the case of the Euclidean norm convergence rate of ``\mathcal{O}\left(k^{2 /(1-m)}\right)``.

"""

# ╔═╡ cf460cb7-9367-4b7b-b2bb-fdb76c1e70f6
md"## Unbounded Problems"

# ╔═╡ 1b20ac49-8453-4f55-89d8-b1e50b3cd85d
cm""""
- As far as I know, for some __unbounded problems__ such __polyhedral approximations__ do not exist.
"""

# ╔═╡ e4ed7740-d159-4d35-a50a-038cd7ac184b
# let
# 	# max(tr_corner_x,max_x)
# 	xc,yc, r = 2.0,2.0,1.0
# 	tr_corner_x,tr_corner_y=(6.4,5.0)
# 	# X = [1.0  2.0  1.29289  1.6398   1.06712  1.55167
#  # 2.0  1.0  1.29289  1.06712  1.6398   1.10613]
# 	d = [1.0;1.0]
# 	θ_t = π/4
# 	v_x = xc + r*cos(θ_t+π)
# 	v_y = yc + r*sin(θ_t+π)
# 	v = [v_x;v_y]
# 	b = dot(d,v)
# 	xs = [1;(b .- d[2]*1) / d[1]]
# 	ys = [(b .- d[1]*1) / d[2];1]
# 	X = [xs';ys']
# 	p = sortperm(X[1, :])  
# 	X_sorted=X[:,p]
# 	# p1 = draw_outer([v[1]],[v[2]],tr_corner=(6.4,5.0))
# 	p1 = draw_outer(X_sorted[1,:],X_sorted[2,:],tr_corner=(6.4,5.0))
# 	p = draw_example1(xc,yc, r; 
# 					  upper_image=false,
# 					  show_center=false,
# 					  tr_corner=(tr_corner_x,tr_corner_y),
# 					  initial_plot=p1
# 					 )
	
# 	# 5) add the two text labels
# 	annotate!(p, tr_corner_x-0.5, tr_corner_y-0.5, text(L"P_1", :center, 12))
	
# end

# ╔═╡ a2d1d14d-c984-4cc6-9fc2-76365f3d9c2f
# let
# 	# max(tr_corner_x,max_x)
# 	xc,yc, r = 2.0,2.0,1.0
# 	tr_corner_x,tr_corner_y=(6.4,5.0)
# 	X = [1.0  2.0  1.29289  1.6398   1.06712  1.55167
#  2.0  1.0  1.29289  1.06712  1.6398   1.10613]
# 	p = sortperm(X[1, :])  
# 	X_sorted=X[:,p]
# 	p1 = draw_outer(X_sorted[1,:],X_sorted[2,:],tr_corner=(6.4,5.0))
# 	p = draw_example1(xc,yc, r; 
# 					  upper_image=false,
# 					  show_center=false,
# 					  tr_corner=(tr_corner_x,tr_corner_y),
# 					  initial_plot=p1
# 					 )
	
# 	# 5) add the two text labels
# 	annotate!(p, tr_corner_x-0.5, tr_corner_y-0.5, text(L"P_5", :center, 12))
	
# end

# ╔═╡ a631fd0c-b9b6-444a-8af9-7dbe8bc7fda4
let
	direction_selection = [
	    :FixedDirection,
	    :FixedPointDirection,
	    :IdealPointDirection,
	    :ApproximateAdjDirection,
	  ]
 	vertex_selection = 
		[
			:RandomVertex, 
			:IdealPointVertex, 
			:InnerPointVertex
		]
	# examp1_option = ExampleOptions(2, 2, 0.05, 100, 5, 50, 30)
	# dims = [2]
	# res1 = example1(
	#   dims,
	#   examp1_option;
	#   direction_selection = [:FixedDirection],
	#   vertex_selection = [:InnerPointVertex],
	# )
	# res1[1].X
	X = [1.0  2.0  1.29289  1.6398   1.06712  1.55167
 2.0  1.0  1.29289  1.06712  1.6398   1.10613]
	md""
end

# ╔═╡ 34c3882f-2178-40f8-80b8-e3f019e5284e
begin
	# Hex-based palette
	dark_palette = Dict(
    :midnight_blue     => "#191970",
    :dark_slate_gray   => "#2F4F4F",
    :dark_olive_green  => "#556B2F",
    :dark_green        => "#006400",
    :dark_slate_blue   => "#483D8B",
    :dark_red          => "#8B0000",
    :dark_magenta      => "#8B008B",
    :charcoal          => "#2E2E3A"
)
	md""
end

# ╔═╡ d71e0ae6-a621-4738-9683-443b4043e371
begin
	struct ExampleOptions
	  pdim::Int
	  ddim::Int
	  ϵ::Float64
	  maxiters::Int
	  trials::Int
	  maxcard::Union{Nothing,Int}
	  time_limit::Union{Nothing,Int}
	end
	ExampleOptions(pdim::Int) = ExampleOptions(pdim, pdim, 0.05, 1000, 1, nothing, nothing)
	ExampleOptions(pdim::Int, ddim::Int) = ExampleOptions(pdim, ddim, 0.05, 1000, 1, nothing, nothing)
	ExampleOptions(eo::ExampleOptions, pdim::Int, ddim::Int) =
	  ExampleOptions(pdim, ddim, eo.ϵ, eo.maxiters, eo.trials, eo.maxcard, eo.time_limit)
end

# ╔═╡ 3c5c8346-86f7-4309-a476-43a35a50f8bc
begin
	

function example1(
  dims::Vector{Int},
  optionts::ExampleOptions;
  direction_selection::Vector{Symbol} = [
    :FixedDirection,
    :FixedPointDirection,
    :IdealPointDirection,
    :ApproximateAdjDirection,
  ],
  vertex_selection::Vector{Symbol} = [:RandomVertex, :IdealPointVertex, :InnerPointVertex],
)
  results = map(dims) do d
    eo = ExampleOptions(optionts, d, d)
    example1(eo, direction_selection = direction_selection, vertex_selection = vertex_selection)
  end
  vcat(results...)
end

function example1(
  optionts::ExampleOptions;
  direction_selection::Vector{Symbol} = [
    :FixedDirection,
    :FixedPointDirection,
    :IdealPointDirection,
    :ApproximateAdjDirection,
  ],
  vertex_selection::Vector{Symbol} = [:RandomVertex, :IdealPointVertex, :InnerPointVertex],
)
  pdim, ddim, ϵ, maxiters, trials, time_limit, maxcard = optionts.pdim,
  optionts.ddim,
  optionts.ϵ,
  optionts.maxiters,
  optionts.trials,
  optionts.time_limit,
  optionts.maxcard
  direction_selection = symbolToDirection.(direction_selection)
  vertex_selection = symbolToVertex.(vertex_selection)
  solver = optimizer_with_attributes(
    Ipopt.Optimizer,
    MOI.Silent() => true,
    "sb" => "yes",
    "max_iter" => 10_000,
  )
  # (1) The function
  f(x) = x
  # (2) The matrix whose columns generate C
  Z = Float64.(Matrix(I(ddim)))  # generating C
  C = NMCone(Z)
  # (3) The matrix whose columns generate C+
  Zplus = copy(Z) # generating C+
  Cp = NMCone(Z)

  # (4) direction inside C
  d = ones(ddim)
  d = d / norm(d)

  # (5) Pascoletti-Serafini scalarization
  function ps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, z)
    @variable(model, x[1:pdim] >= 0)
    @NLconstraint(model, sum((x[i] - 2)^2 for i = 1:pdim) <= 1)
    @constraint(model, Zplus * (v + z * d - f(x)) .>= 0)
    @objective(model, Min, z)
    optimize!(model)
    if (has_values(model))
      xv = value.(x)
      return xv, value.(z)
    else
      return nothing
    end
  end

  # (6) To find the realized approximation error is dH(P K, P)
  function ps2(v::Vector{<:Number})
    model = Model(solver)
    @variable(model, z[1:ddim])
    @variable(model, x[1:pdim] >= 0)
    @NLconstraint(model, sum((x[i] - 2)^2 for i = 1:pdim) <= 1)
    @constraint(model, Zplus * (v + z - f(x)) .>= 0)
    @NLobjective(model, Min, sum(z[i]^2 for i = 1:ddim))
    optimize!(model)
    if (has_values(model))
      return norm(value.(z))
    else
      return nothing
    end
  end

  # (7) The Lagrange dual of the PS problem
  function dps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, t)
    @variable(model, x[1:pdim] >= 0)
    @variable(model, w[1:ddim])
    @constraint(model, dot(w, d) == 1)
    @constraint(model, Z * w .>= 0)
    @NLconstraint(model, sum((x[i] - 2)^2 for i = 1:pdim) <= 1)
    @constraint(model, dot(f(x), w) >= t)
    @objective(model, Max, t - dot(w, v))
    optimize!(model)
    if (has_values(model))
      return value.(w)
    else
      return nothing
    end
  end

  # (8) the weighted sum scalarization
  function ws(w::Vector{<:Number})
    model = Model(solver)
    @variable(model, x[1:pdim] >= 0)
    @NLconstraint(model, sum((x[i] - 2)^2 for i = 1:pdim) <= 1)
    @objective(model, Min, dot(w, f(x)))
    optimize!(model)
    if (has_values(model))
      return value.(x)
    else
      return nothing
    end
  end

  # (9) Get the ideal point using the weighted sum scalarization
  function getIdealPoint()
    yI = map(1:ddim) do j
      model = Model(solver)
      @variable(model, x[1:pdim] >= 0)
      @NLconstraint(model, sum((x[i] - 2)^2 for i = 1:pdim) <= 1)
      @objective(model, Min, f(x)[j])
      optimize!(model)
      if (has_values(model))
        xv = value.(x)
        return xv[j]
      else
        return nothing
      end
    end
    return yI
  end

  # (10) The ideal point   
  yI = abs.(getIdealPoint())
  # (11) The problem 
  problem = MOAProblem(pdim, ddim, C, Cp, f, ws, ps, dps, d, yI)

  # array to hold the results
  results = Vector{MOASuccessResult}(undef, length(direction_selection) * length(vertex_selection))

  i = 0
  for vertex in vertex_selection
    tries, average_it = (Symbol(vertex) == :RandomVertex) ? (trials, true) : (1, false)
    printstyled("Staring $(String(Symbol(vertex))) with $tries tries.\n", color = :light_magenta)
    printstyled("=========================================\n", color = :blue)
    tries
    for direction in direction_selection
      printstyled(
        "---- with $(String(Symbol(direction))) direction selection .\n",
        color = :light_magenta,
      )
      total_trial_counter = 0
      success_trial_counter = 0
      Xs = Vector{MOASuccessResult}(undef, tries)
      while true
        if total_trial_counter > 20 * tries
          printstyled("Could not find $tries successful results.\n", color = :red)
          break
        end
        if success_trial_counter == tries
          break
        end
        t1 = Dates.now()
        options = MOAOptions(
          ϵ,
          maxiters,
          (c, card) -> begin
            cond_iters = c >= maxiters
            cond_card = isnothing(maxcard) ? false : card >= maxcard
            cond_time = if isnothing(time_limit)
              false
            else
              t2 = Dates.now()
              return round(t2 - t1, Second) > Second(time_limit)
            end
            return cond_iters || cond_card || cond_time
          end,
        )
        X, t = @timed paa(direction, vertex, problem, options)
        if isa(X, MOASuccessResult)
          success_trial_counter += 1
          printstyled(
            "Got a successfull result $success_trial_counter with message $(X.message).\n",
            color = :blue,
          )
          V, = getVertices(X.Pk)
          hdist1 = hausdorffDistance(V, ps2)
          # hdist2 = hausdorffDistance(X.Vk, PS_distance(Z, f, n, n, solver))
          # printstyled("distance 1 = $hdist1 and distance 2 = $hdist2\n")
          X = MOASuccessResult(X, String(Symbol(direction)), String(Symbol(vertex)), t)
          Xs[success_trial_counter] = MOASuccessResult(X, t, hdist1)

        end
        total_trial_counter += 1
      end
      if average_it
        name = MOAResultName(String(Symbol(direction)), String(Symbol(vertex)))
        R = MOASuccessResult(
          name,
          Xs[tries].message,
          Xs[tries].X,
          Xs[tries].Vk,
          Xs[tries].Pk,
          mean([r.SC for r in Xs]),
          mean([r.Card for r in Xs]),
          mean([r.T for r in Xs]),
          mean([r.HD for r in Xs]),
          problem,
        )
        i = i + 1
        results[i] = R
      else
        results[i+1:i+tries] = Xs
        i = i + tries
      end
    end
  end

  return results
end
end

# ╔═╡ d899a2e5-ed49-4bee-9dad-a7934ea5d560
begin
	function draw_example1(xc=2.0, yc=2.0, radius=1.0; 
						   upper_image=true, 
						   tr_corner=(5.0,5.0),
						   show_center=true,
						   pareto_front=false,
						   initial_plot=plot(),
						   title=""
						  )
			# Pick a theme
		theme(:wong2)
		tr_corner_x,tr_corner_y = tr_corner
		
		# Parameterize the boundary circle
		θ = LinRange(0, 2π, 200)
		x = xc .+ radius .* cos.(θ)
		y = yc .+ radius .* sin.(θ)
		p = if upper_image
			rect_x1 = [xc, tr_corner_x, tr_corner_x, xc]
			rect_y1 = [yc-radius, yc-radius, tr_corner_y, tr_corner_y]
			rect_x2 = [xc-radius, tr_corner_x, tr_corner_x, xc-radius]
			rect_y2 = [yc, yc, tr_corner_y, tr_corner_y]
			the_cone = [
				Shape(rect_x1,rect_y1),
				Shape(rect_x2,rect_y2),
			]
			plot(the_cone,
					 seriestype   = :shape,
					 fillcolor    =:darkgrey, 
					 fillalpha    = 1.0,
					 linecolor    = :transparent,
					)
			else
				initial_plot
			end
		# Draw the filled disk
		p = plot(p,
		    x, y;
		    seriestype   = :shape,      # draw as a filled polygon
		    fillcolor    = :gray,    # fill color
		    fillalpha    = 1.0,         # transparency
		    linecolor    = :gray,       # boundary color
		    linewidth    = 1.5,
			aspect_ratio = :equal,      # equal scaling on x and y
		    legend       = false,
		    xlabel       = "x",
		    ylabel       = "y",
		    title        = title,
			frame_style  = :origin
		)
		if pareto_front
		θ = LinRange(π, 3π/2, 100)
		x = xc .+ radius .* cos.(θ)
		y = yc .+ radius .* sin.(θ)
		scatter!(
		    x, y;
		    color    	 = :green,       # boundary color
		    markersize 	 = 2,
		    aspect_ratio = :equal,      # equal scaling on x and y
		    legend       = false,
		    # xlabel       = none,
		    # ylabel       = nothing,
		    frame_style  = :origin
		)
		end
		# Optionally mark the center
		if show_center
		scatter!([xc], [yc];
		    color   = :red,
		    marker  = :circle,
		    markersize = 2,
		)
		annotate!(p, xc, yc, text(L"(2,2)", :center, 14))
		end
		p
	end
	function draw_outer(xs,ys;tr_corner=(5.0,5.0))
		min_x, max_x = minimum(xs), maximum(xs)
		min_y, max_y = minimum(ys), maximum(ys)
		tr_corner_x,tr_corner_y = tr_corner
		Pk = Shape(
			[xs...,max(max_x,tr_corner_x),max(max_x,tr_corner_x),min_x],
			[ys...,min_y,tr_corner_y,tr_corner_y]
		)
		p = plot(Pk,
				seriestype   = :shape,
				fillcolor    =:darkgrey, 
				fillalpha    = 1.0,
				linecolor    = :transparent,
				label = ""
		)
		p = scatter(p,xs,ys,
				 label = ""
				   )
		p
	end
end

# ╔═╡ 9f14a6ae-8756-4059-a984-3b935ff598a0
let
	draw_example1(upper_image=false,title=L"f(\Omega)=\Omega",pareto_front=true)
end	

# ╔═╡ d9c32ddb-9922-4b73-882d-bc353e0dd761
let
	xc,yc, r = 2.0,2.0,1.0
	p = draw_example1(xc,yc, r;show_center=false,title="Upper Image")
	annotate!(p, xc, yc+0.5, text(L"f(\mathcal{\Omega})", :center, 14))
	annotate!(p, xc+0.6, yc + 1.2, text(L"\mathcal{P}=f(\Omega)+\mathbb{R}^2_+",    :left,   14))
			
end

# ╔═╡ a707080d-33ef-4dda-b2a3-72a23dc90db3
let
	# pick a theme
	# theme(:wong2)
	xc,yc, r = 2.0,2.0,1.0
	tr_corner_x,tr_corner_y=(6.4,5.0)
	p = draw_example1(xc,yc, r; 
					  show_center=false,
					  tr_corner=(tr_corner_x,tr_corner_y)
					 )
	
	# 5) add the two text labels
	annotate!(p, xc, yc+0.5, text(L"f(\mathcal{\Omega})", :center, 12))
	
	# 6) inner P
	xs = [tr_corner_x, xc, xc-1/sqrt(2), xc-r, xc-r]
	ys = [yc-r, yc-r, xc-1/sqrt(2), yc, tr_corner_y]
	plot!(
		    p, xs, ys,
		    # seriestype   = :shape,
		    # fillcolor    = :blue,
			marker    	 = :circle,
			markersize 	 = 3.5,
			markercolor  = "#483D8B",
		    fillalpha    = 1,
		    linecolor    = "#483D8B",
		    linewidth    = 1.5,
			linestyle    = :dash,
		    label        = "",
		)

	# 7) outer P
	xs = [tr_corner_x, xc, xc-1/sqrt(2), xc-r, xc-r]
	ys = [yc-r, yc-r, xc-1/sqrt(2), yc, tr_corner_y]
	ϵ = 0.1
	plot!(
		    p, xs .- ϵ, ys .- ϵ,
		    # seriestype   = :shape,
		    # fillcolor    = :blue,
			marker    	 = :circle,
			markersize 	 = 1.0,
			markercolor  = dark_palette[:dark_magenta],
		    fillalpha    = 1,
		    linecolor    = dark_palette[:dark_magenta],
		    linewidth    = 1.5,
			linestyle    = :dash,
		    label        = "",
		)
	annotate!(p, xc+0.4, yc+1.8, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+)",    :left,   12,"#483D8B"))	
	annotate!(p, xc+0.4, yc+1.4, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+ -\mathbb{1}\{\epsilon\})",    :left,   12,dark_palette[:dark_magenta]))	
end


# ╔═╡ 5b0e766a-28fd-467d-8c70-717730ec47ef
let
	# pick a theme
	xc,yc, r = 2.0,2.0,1.0
	tr_corner_x,tr_corner_y=(6.4,5.0)
	p = draw_example1(xc,yc, r; 
					  show_center=false,
					  tr_corner=(tr_corner_x,tr_corner_y)
					 )
	
	# 5) add the two text labels
	annotate!(p, xc, yc+0.5, text(L"f(\mathcal{\Omega})", :center, 12))
	
	# 6) inner P
	xs = [tr_corner_x, xc, xc-1/sqrt(2), xc-r, xc-r]
	ys = [yc-r, yc-r, xc-1/sqrt(2), yc, tr_corner_y]
	plot!(
		    p, xs, ys,
		    # seriestype   = :shape,
		    # fillcolor    = :blue,
			marker    	 = :circle,
			markersize 	 = 1.0,
			markercolor  = "#483D8B",
		    fillalpha    = 1,
		    linecolor    = "#483D8B",
		    linewidth    = 1.5,
			linestyle    = :dash,
		    label        = "",
		)
	scatter!(
		    p, xs[2:4], ys[2:4],
		    # seriestype   = :shape,
		    # fillcolor    = :blue,
			marker    	 = :circle,
			markersize 	 = 6.0,
			markeralpha  = 0.8,
			markercolor  = "#483D8B",
		    # fillalpha    = 0.2,
		    # linecolor    = "#483D8B",
		    label        = "",
		)
	# 7) outer P
	xs = [tr_corner_x, xc, xc-1/sqrt(2), xc-r, xc-r]
	ys = [yc-r, yc-r, xc-1/sqrt(2), yc, tr_corner_y]
	ϵ = 0.1
	plot!(
		    p, xs .- ϵ, ys .- ϵ,
		    # seriestype   = :shape,
		    # fillcolor    = :blue,
			marker    	 = :circle,
			markersize 	 = 1.0,
			markercolor  = dark_palette[:dark_magenta],
		    fillalpha    = 1,
		    linecolor    = dark_palette[:dark_magenta],
		    linewidth    = 1.5,
			linestyle    = :dash,
		    label        = "",
		)
	annotate!(p, xc+0.4, yc+1.8, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+)",    :left,   12,"#483D8B"))	
	annotate!(p, xc+0.4, yc+1.4, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+ +\mathbf{B}(\mathbf{0},\epsilon))",    :left,   12,dark_palette[:dark_magenta]))
	
	
end


# ╔═╡ a5bde7fe-472a-46e5-9cab-494e9c945493
let
	xc,yc, r = 2.0,2.0,1.0
	tr_corner_x,tr_corner_y=(6.4,5.0)
	p = draw_example1(xc,yc, r; 
					  show_center=false,
					  tr_corner=(tr_corner_x,tr_corner_y)
					 )
	
	annotate!(p, xc, yc+0.5, text(L"f(\mathcal{\Omega})", :center, 12))
	θ_t = π/4
	p_x = xc + r*cos(θ_t+π)
	p_y = yc + r*sin(θ_t+π)
	d = [cos(θ_t), sin(θ_t)]  # points toward the center
	quiver!(p, [0.0], [0.0],
	    quiver      = ([d[1]], [d[2]]),
	    arrow       = true,
	    linewidth   = 1.2,
	    linecolor   = :blue,
	    linestyle   = :dashdot,
	)
	annotate!(p, p_x-0.5,  p_y-0.5,  text(L"d",             :blue, 14))

	scatter!(p, [0.0], [0.0];
	    color      = :black,
	    markersize = 4,
	)
	annotate!(p, -0.4,  -.2,  text(L"v=\mathbf{0}",             14))
	b = dot(d, [p_x, p_y])
	xs = [0.0;b/d[1]]
	ys = (b .- d[1]*xs) ./ d[2]
	plot!(p, xs, ys;
	    linecolor  = :red,
	    linestyle  = :dash,
	    linewidth  = 2,
	)
	
	pt_x, pt_y = b/d[1], b/d[2]
	scatter!(p, [pt_x, 0.0], [0.0, pt_y];
	    color      = :red,
	    markersize = 4,
	)
	
	annotate!(p, xc+2.8, yc+2.8, text(L"\mathcal{P}",    :left,   12))
	annotate!(p, -0.3,  pt_y,  text(L"H",      :red,       14))
end

# ╔═╡ 7a491a11-2b36-41f5-84d7-ed02d8d5326e
let
	# max(tr_corner_x,max_x)
	
	xc,yc, r = 2.0,2.0,1.0
	tr_corner_x,tr_corner_y=(6.4,5.0)
	X,Ptitle = if steps=="3"
	X1 = [1.0  2.0  1.29289  1.6398   1.06712  1.55167
 2.0  1.0  1.29289  1.06712  1.6398   1.10613]
		X1,L"P_5"
	elseif steps == "2"	
		d = [1.0;1.0]
	θ_t = π/4
	v_x = xc + r*cos(θ_t+π)
	v_y = yc + r*sin(θ_t+π)
	v = [v_x;v_y]
	b = dot(d,v)
	xs = [1;(b .- d[2]*1) / d[1]]
	ys = [(b .- d[1]*1) / d[2];1]
	X1 = [xs';ys']
		X1,L"P_1"
	else
		[1.0;1.0],L"P_0"
	end
	steps
	p = sortperm(X[1, :])  
	X_sorted=X[:,p]
	p1 = draw_outer(X_sorted[1,:],X_sorted[2,:],tr_corner=(6.4,5.0))
	p = draw_example1(xc,yc, r; 
					  upper_image=false,
					  show_center=false,
					  tr_corner=(tr_corner_x,tr_corner_y),
					  initial_plot=p1
					 )
	
	# 5) add the two text labels
	annotate!(p, tr_corner_x-0.5, tr_corner_y-0.5, text(L"%$(Ptitle)", :center, 12))
	
end

# ╔═╡ 66ffc650-93e6-4705-8212-7958f75bb1e1
begin
    function add_space(n=1)
        repeat("&nbsp;", n)
    end
    function post_img(img::String, w=500)
        res = Resource(img, :width => w)
        cm"""
      <div class="img-container">

      $(res)

      </div>"""
    end
    function poolcode()
        cm"""
      <div class="img-container">

      $(Resource("https://www.dropbox.com/s/cat9ots4ausfzyc/qrcode_itempool.com_kfupm.png?raw=1",:width=>300))

      </div>"""
    end
    function define(t="")
        beginBlock("Definition", t)
    end
    function remark(t="")
        beginBlock("Remark", t)
    end
    function remarks(t="")
        beginBlock("Remarks", t)
    end
    function bbl(t)
        beginBlock(t, "")
    end
    function bbl(t, s)
        beginBlock(t, s)
    end
    ebl() = endBlock()
    function theorem(s)
        bth(s)
    end
    function bth(s)
        beginTheorem(s)
    end
    eth() = endTheorem()
    ex(n::Int; s::String="") = ex("Example $n", s)
    ex(t::Int, s::String) = example("Example $t", s)
    ex(t, s) = example(t, s)
    function beginBlock(title, subtitle)
        """<div style="box-sizing: border-box;">
       	<div style="display: flex;flex-direction: column;border: 6px solid rgba(200,200,200,0.5);box-sizing: border-box;">
       	<div style="display: flex;">
       	<div style="background-color: #FF9733;
       	    border-left: 10px solid #df7300;
       	    padding: 5px 10px;
       	    color: #fff!important;
       	    clear: left;
       	    margin-left: 0;font-size: 112%;
       	    line-height: 1.3;
       	    font-weight: 600;">$title</div>  <div style="olor: #000!important;
       	    margin: 0 0 20px 25px;
       	    float: none;
       	    clear: none;
       	    padding: 5px 0 0 0;
       	    margin: 0 0 0 20px;
       	    background-color: transparent;
       	    border: 0;
       	    overflow: hidden;
       	    min-width: 100px;font-weight: 600;
       	    line-height: 1.5;">$subtitle</div>
       	</div>
       	<p style="padding:5px;">
       """
    end
    function beginTheorem(subtitle)
        beginBlock("Theorem", subtitle)
    end
    function endBlock()
        """</p></div></div>"""
    end
    function endTheorem()
        endBlock()
    end
    ex() = example("Example", "")
    function example(lable, desc)
        """<div class="example-box">
    <div class="example-header">
      $lable
    </div>
    <div class="example-title">
      $desc
    </div>
    <div class="example-content">
      
  </div>
		"""
    end
	md"``\color{white}{\text{Packages}}``"
end

# ╔═╡ f23b7d03-2a80-492f-b933-317f735ec553
cm"""
$(bbl("Problem","Convex Vector Optimization Problem"))

We consider the following optimization problem
```math
\begin{array}{ll}
\min_C & f(x) \\
\textrm{subject to}
&  x \in \Omega
\end{array}\tag{P}
```

where
- ``\displaystyle f:\mathbb{R}^n \to \mathbb{R}^m`` is __convex__.
- ``\displaystyle \Omega \subseteq \mathbb{R}^n`` is __convex__.
- ``\displaystyle C \subseteq \mathbb{R}^m`` is a __closed, convex, pointed__ cone with __nonempty interior__.

Therefore this problem is called __Convex Vector Optimization Problem__.

Furthermore, "``\min_C``"" in the sense that if ``x, y \in \mathbb{R}^n``, then
```math
f(x) \leq_C f(y) \quad \Rightarrow \quad f(y)-f(x) \in C.
```



"""

# ╔═╡ 57b118d7-a7dd-4f57-95e3-ef7d3fc334bf
cm"""
$(ex("","Example"))
For demonstration, we consider the following example. Let ``x\in\mathbb{R}^2`` and solve 
```math
\begin{array}{ll}
\min_{\mathbb{R}^2_+} & x \\
\textrm{subject to}
&  \|x-\mathbb{2}\|_2^2 \leq 1
\end{array}
```

"""

# ╔═╡ 8a7307cb-0ad4-4e46-b18d-92eb89ac82cb
cm"""
$(define("C-Convex")) 
A function ``f: \mathbb{R}^n \rightarrow \mathbb{R}^p`` is said to be __C-convex__ if 

```math
f(\lambda x+(1-\lambda) y) \leq_C \lambda f(x)+(1-\lambda) f(y)
``` 
for all ``x, y \in \mathbb{R}^n, \lambda \in[0,1]``.
"""

# ╔═╡ 871cf958-6093-4004-99fd-ead3a2561e50
cm"""
$(define("Minimal Elements")) 
For a set ``S \subseteq \mathbb{R}^m``, 
- a point ``s \in S`` is a __``C``-minimal element of ``S``__ if ``(\{s\}-C \backslash\{0\}) \cap S=\emptyset``. 
- If ``C`` has nonempty interior and ``(\{s\}-\operatorname{int} C) \cap S=`` ``\emptyset``, then ``s \in S`` is called a __weakly ``C``-minimal element of ``S``__. 

The set of all ``C``-minimal (weakly ``C``-minimal) elements of ``S`` is denoted by ``\operatorname{Min}_C S\left(\mathrm{wMin}_C S\right)``.
"""

# ╔═╡ 647e1f03-9a12-4239-b649-8d9aabbaa7c2
cm"""
$(bbl("Assumption",""))

```math
C \quad \text{is a polyhedral}.
```
So its dual cone is given by 
```math
C^{+}:=\operatorname{co conv}\left\{w^1, \ldots, w^l\right\}, 
``` 
and ``\mathcal{\Omega}`` is compact with nonempty interior.
$(ebl())

$(define("Polyhedral Convex Set"))
A set ``S\in\mathbb{R}^m`` is a __polyhedral convex set__ if it is of the form 
```math
S=\left\{y \in \mathbb{R}^m \mid A^T y \geq b\right\},\quad \color{red}{\text{Halfspace representation (H-representation)}}
```
where ``A \in`` ``\mathbb{R}^{m \times n}, b \in \mathbb{R}^n``. 
$(ebl())

__Note__:

``S`` can also be represented as 
```math
S = \operatorname{conv} V + \operatorname{co} D\quad \color{red}{\text{Vertix representation (V-representation)}}
```
where ``V \subseteq \mathbb{R}^m`` and ``D \subseteq \mathbb{R}^m`` are the finite sets of vertices and extreme directions of ``S``, respectively.

$(define("Upper Image"))
The upper image for probem (P) is defined as 
```math
\mathcal{P}:=\operatorname{cl}(f(\mathcal{\Omega})+C)
```
$(ebl())


"""

# ╔═╡ 3bf3f0c2-496f-41b5-b2bc-0c0cb7b6feab
cm"""

$(ex("Remarks",""))
- Note that 
```math
\begin{array}{lcl}
\mathbb{R}^m_+&=&\left\{y \in \mathbb{R}^m \mid I y \geq \mathbf{0}\right\},\\
&=&\left\{\mathbf{0}\right\}+\operatorname{co}\left\{e_1,e_2,\cdots,e_m\right\},\\
\end{array}
```
- ``\mathcal{P}`` is a closed convex set and under the assumptions of the problem, 
```math
\mathcal{P}=f(\mathcal{\Omega})+C.
```
- the set of all weakly ``C``-minimal points of ``\mathcal{P}`` is ``\operatorname{bd}\mathcal{P}``.

- Since ``\Omega`` is compact, we know that problem (P) is __bounded__ in the sense that 
```math
\mathcal{P} \subseteq \{y\} + C, \quad \text{for some}\quad y \in \mathbb{R}^m.
```

"""

# ╔═╡ bb397f1a-1cde-4e83-afe3-e9dbad4a485f
cm"""
$(define("Dual Cone"))
Let ``C \subset \mathbb{R}^m`` be a nonempty closed convex cone. The __dual cone__ of ``C`` is defined as 
```math
C^{+}:=\left\{a \in \mathbb{R}^m \mid \forall x \in C: a^{\top} x \geq 0\right\}
``` 
$(ebl())

$(bbl("Remarks",""))
- ``C^{+}`` is a closed convex cone. 
- If ``C \subseteq \mathbb{R}^m`` is a polyhedral cone, then ``C^{+}``is also polyhedral and can be written as 

```math
C^{+}=\operatorname{co conv}\left\{w^1, \ldots, w^l\right\},
``` 
$(add_space(10))where ``w^1, \ldots, w^l \in \mathbb{R}^m`` are the extreme directions of ``C^{+}``for some ``l \geq 0``. 

$(add_space(10))In this case, 
```math
y^1 \leq_C y^2\; \text{holds if and only if } \left(w^i\right)^{\top} y^1 \leq\left(w^i\right)^{\top} y^2 \quad \text{for all } i \in\{1, \ldots, l\}.
```
"""

# ╔═╡ 035fe90d-6c65-4a6a-b345-10d44e12b131
cm"""
$(define("ϵ-solution"))
Let ``c \in \operatorname{int} C`` be fixed. For ``\epsilon>0``, a nonempty finite set ``\overline{\mathcal{\Omega}} \subseteq \mathcal{\Omega}`` of (weak) minimizers is a finite __(weak) ``\epsilon``-solution__ with respect to ``c`` if 
```math
\text{conv }f(\overline{\mathcal{\Omega}})+ C-\epsilon\{c\} \supseteq \mathcal{P}.
```
$(ebl())

"""

# ╔═╡ 58d3a13b-87cd-481c-a208-7259607b034c
cm"""
$(define("ϵ-solution, no direction")) 
For ``\epsilon>0``, a nonempty finite set ``\overline{\mathcal{\Omega}} \subseteq \mathcal{\Omega}`` of (weak) minimizers is a finite (weak) ``\epsilon``-solution if ``\operatorname{conv} f(\overline{\mathcal{\Omega}})+C+B(0, \epsilon) \supseteq \mathcal{P}``.
"""

# ╔═╡ 6c7366ee-895d-46b6-ba30-7604ba6213a9
cm"""
$(ex("Methods of Solutions",""))

- There are different solution approaches to solve __(P)__. 

- The __main idea__ of these approaches is to __generate (weak) minimizers__ for (P) in a structured way. 

- One way of generating (weak) minimizers is to __solve scalarization models__

"""

# ╔═╡ f72ebc57-eb4e-457a-926f-27317caffbeb
cm"""
$(bbl("Proposition",""))
An optimal solution ``x \in \mathcal{X}`` of ``(\mathrm{WS}w)`` is a __weak minimizer of ``(\mathrm{P})``__ if ``w \in C^{+} \backslash\{0\}``. 

__Conversely__, for any weak minimizer ``x \in \mathcal{X}``, there exists ``w \in C^{+} \backslash\{0\}`` such that ``x`` is an optimal solution of ``(\mathrm{WS}w)``.
"""

# ╔═╡ 7b3c1689-867e-47eb-9641-107d5a866348
cm"""
$(bbl("Proposition",""))
If ``\left(x, z\right) \in \mathbb{R}^{n+1}`` is an optimal solution of problem ``(\operatorname{PS}(v, d))``, then ``x`` is a weak minimizer. Moreover, ``v+z d \in \operatorname{bd} \mathcal{P}``.
$(ebl())

$(ex("Remark",""))
- Note that ``f(x) \leq_C v+z d`` holds for some ``x \in \mathcal{\Omega}`` if and only if ``v+z d \in \mathcal{P}``. 

To see, assume ``f(x) \leq_C v+z d`` holds for some ``x \in \mathcal{\Omega}``, that is, ``v+z d-f(x) \in C`` holds. Then, we have ``v+z d \in\{f(x)\}+C \subseteq f(\mathcal{\Omega})+C=\mathcal{P}``. 

The other implication follows similarly.
"""

# ╔═╡ cca4e694-eca3-46d3-ac0b-19392f9060a6
cm"""
$(bbl("Proposition","Existence of solutions to PS(v,d) and DPS(v,d)"))
Let ``v \in \mathbb{R}^m``. If ``d \in \operatorname{int} C``, then there exist optimal solutions to ``(\operatorname{PS}(v, d))`` and (DPS(v,d)). Moreover, the optimal values of the two problems __coincide__.
"""

# ╔═╡ e63032ba-c53b-4e92-91eb-f611ceb89903
cm"""
$(bbl("Proposition",""))
For ``v \notin \mathcal{P}, y \in \operatorname{int} \mathcal{P}`` and ``d=y-v``, both ``(\operatorname{PS}(v, d))`` and (DPS ``(v, d))`` have optimal solutions and the optimal values coincide.
"""

# ╔═╡ 39da7d26-392c-4d69-b4c9-59876257c657
cm"""

$(bbl("Proposition","Existence of Supporting hyperplane to the upper image"))
Let ``v, d \in \mathbb{R}^m`` and ``\left(x^*, z^*\right)``, ``w^* \in \mathbb{R}^m`` be the optimal solutions for ``(\operatorname{PS}(v, d))`` and its Lagrange dual, respectively. 

If the optimal objective function values of ``(\operatorname{PS}(v, d))`` and ``(\mathrm{DPS}(v, d))`` coincides, then 
```math
H:=\left\{y \in \mathbb{R}^m \mid\left(w^*\right)^{\top} y=\right.
\left.\left(w^*\right)^{\top} v+z^*\right\}
```
is a __supporting hyperplane__ for ``\mathcal{P}`` at ``y^*=v+z^* d`` and 

```math
\mathcal{H}=\left\{y \in \mathbb{R}^m \mid\right. \left.\left(w^*\right)^{\top} y \geq\left(w^*\right)^{\top} v+z^*\right\} \supseteq \mathcal{P}
``` 
"""

# ╔═╡ f8ee59ea-0832-4d0e-bd51-50bf8b63e222
cm"""
$(ex("Remarks","Some Geometry"))
Let ``S \subseteq \mathbb{R}^m`` be a convex set. 

> A __hyperplane__ given by ``\left\{y \in \mathbb{R}^m \mid a^{\top} y=b\right\}`` for some ``a \in \mathbb{R}^m \backslash\{0\}, b \in \mathbb{R}`` is a __supporting hyperplane__ of ``S`` if ``S \subseteq\left\{y \in \mathbb{R}^m \mid a^{\top} y \geq b\right\}`` and there exists ``s \in S`` with ``a^{\top} s=b``. 

> A convex subset ``F \subseteq S`` is called a __face__ of ``S`` if ``\lambda x+(1-\lambda) y \in F \text{ with } x, y \in S \text{ and } 0<\lambda<1\Rightarrow x, y \in F.``

> A zero-dimensional face is an __extreme point (or vertex)__.

> A one-dimensional face is an __edge__ of ``S``. 

> A recession direction ``z \in \mathbb{R}^m \backslash\{0\}`` of convex set ``S`` is said to be an __extreme direction__ of ``S`` if ``\left\{v+r z \in \mathbb{R}^m \mid r \geq 0\right\}`` is a face for some extreme point ``v`` of ``S``


> A vector ``z \in \mathbb{R}^m \backslash\{0\}`` is a __recession direction__ of ``S``, if ``y+\gamma z \in S`` for all ``\gamma \geq 0, y \in S``. The set of all recession directions of ``S`` is the recession cone of ``S``.

> The __problem of finding the set of all vertices and extreme directions__ of a polyhedral convex set S, given its halfspace representation is called the __vertex enumeration problem__. In our implementation, we used __`Polyhedra.jl`__ which is julia package for polyhedra manipulations.

"""


# ╔═╡ 7950a01e-48bb-4094-9842-a7dd828ceda6
cm"""
$(bth(""))
When terminates, Previous Algorithm returns a __finite weak ``\epsilon``-solution__.
"""

# ╔═╡ e7657c9d-1ff1-46c8-a19e-99329474cae8
@htl("""
<style>
@import url("https://mmogib.github.io/math102/custom.css");

ul {
  list-style: none;
}

ul li:before {
  content: '▶ ';
}

ul li li:before {
  content: '◼';
}

.p40 {
	padding-left: 40px;
}

example-box {
      max-width: 600px;           /* Limits the box width */
      margin: 2rem auto;          /* Centers the box and adds vertical spacing */
      border: 1px solid #ccc;     /* Light border */
      border-radius: 4px;         /* Slightly rounded corners */
      overflow: hidden;           /* Ensures the box boundary clips its children */
      box-shadow: 0 2px 6px rgba(0, 0, 0, 0.1); /* Subtle shadow */
      font-family: Arial, sans-serif;
    }

    /* Header area for "EXAMPLE 1" */
    .example-header {
      background: linear-gradient(90deg, #cc0000, #990000);
      color: #fff;
      font-weight: bold;
      font-size: 1.1rem;
      padding: 0.75rem 1rem;
      border-bottom: 1px solid #990000;
    }

    /* Sub-header area for the title or subtitle */
    .example-title {
      background-color: #f9f9f9;
      font-weight: 600;
      font-size: 1rem;
      padding: 0.75rem 1rem;
      margin: 0;                  /* Remove default heading margins */
      border-bottom: 1px solid #eee;
    }

    /* Main content area for the mathematical statement or instructions */
    .example-content {
      padding: 1rem;
      line-height: 1.5;
    }

    /* Optional styling for inline math or emphasis */
    em {
      font-style: italic;
      color: #333;
    }
</style>
""")


# ╔═╡ Cell order:
# ╟─1c7f4f9c-0954-408f-927f-66d5929fc006
# ╟─85b8d34e-d89d-49f0-b44e-fc0acee01d73
# ╟─8c4c7c31-e590-4ff5-bc7b-f14d15293980
# ╟─dc6dbd80-444f-4697-ad54-476b16cd5fe6
# ╟─f23b7d03-2a80-492f-b933-317f735ec553
# ╟─c84c7d54-89b8-481e-9b2b-26ede92942f6
# ╟─57b118d7-a7dd-4f57-95e3-ef7d3fc334bf
# ╟─9f14a6ae-8756-4059-a984-3b935ff598a0
# ╟─764259e7-1d76-495b-8553-89f72caf3e46
# ╟─8a7307cb-0ad4-4e46-b18d-92eb89ac82cb
# ╟─871cf958-6093-4004-99fd-ead3a2561e50
# ╟─647e1f03-9a12-4239-b649-8d9aabbaa7c2
# ╟─3bf3f0c2-496f-41b5-b2bc-0c0cb7b6feab
# ╠═d9c32ddb-9922-4b73-882d-bc353e0dd761
# ╟─bb397f1a-1cde-4e83-afe3-e9dbad4a485f
# ╟─035fe90d-6c65-4a6a-b345-10d44e12b131
# ╟─a707080d-33ef-4dda-b2a3-72a23dc90db3
# ╟─58d3a13b-87cd-481c-a208-7259607b034c
# ╠═5b0e766a-28fd-467d-8c70-717730ec47ef
# ╟─1c060d79-1857-495b-994a-90273540afde
# ╟─6c7366ee-895d-46b6-ba30-7604ba6213a9
# ╟─fcb71315-2400-4c64-89ed-a68de52c2cbd
# ╟─d17c513c-768d-4dba-98ae-e7312e84535b
# ╟─f72ebc57-eb4e-457a-926f-27317caffbeb
# ╟─d13e1a2d-bb95-4a45-b920-41abe2ad1a3e
# ╟─204efc5c-9834-4e28-ac2b-26bd1a659a1c
# ╟─7b3c1689-867e-47eb-9641-107d5a866348
# ╟─cca4e694-eca3-46d3-ac0b-19392f9060a6
# ╟─e63032ba-c53b-4e92-91eb-f611ceb89903
# ╟─39da7d26-392c-4d69-b4c9-59876257c657
# ╠═a5bde7fe-472a-46e5-9cab-494e9c945493
# ╟─f9a1dfe1-bc71-42c4-a6e8-3f422eba7507
# ╟─e3f5835c-18fc-4f00-92f5-e55598e68081
# ╟─c6db5f72-a8fa-481c-85e9-7c587e0589b4
# ╟─f8ee59ea-0832-4d0e-bd51-50bf8b63e222
# ╟─574605cb-46c1-4470-b4aa-b4a5e6fe6f8a
# ╟─456f82ae-67de-4050-8d47-e1c106bdd5d1
# ╟─22e078de-0f08-4dd9-bb15-333f5cd7f22b
# ╟─b5d6e08e-318b-4f30-96c3-fe096f6da8dd
# ╟─7a491a11-2b36-41f5-84d7-ed02d8d5326e
# ╟─7950a01e-48bb-4094-9842-a7dd828ceda6
# ╟─453c1441-2b72-4650-a317-56aec0def692
# ╟─f82ea70e-0c04-4275-afbf-4cff107940bf
# ╟─44580f32-df76-484b-9d2f-d09723747e66
# ╟─a9a7ac5b-6f09-4c90-a867-aee556f27c26
# ╟─ef9f9816-74df-4d2d-ac27-12ccd4c1c755
# ╟─513d75b6-63ea-40f4-89b1-f1c79ac83c88
# ╟─cf460cb7-9367-4b7b-b2bb-fdb76c1e70f6
# ╠═1b20ac49-8453-4f55-89d8-b1e50b3cd85d
# ╟─e4ed7740-d159-4d35-a50a-038cd7ac184b
# ╟─a2d1d14d-c984-4cc6-9fc2-76365f3d9c2f
# ╟─a631fd0c-b9b6-444a-8af9-7dbe8bc7fda4
# ╟─34c3882f-2178-40f8-80b8-e3f019e5284e
# ╟─d71e0ae6-a621-4738-9683-443b4043e371
# ╟─3c5c8346-86f7-4309-a476-43a35a50f8bc
# ╠═d899a2e5-ed49-4bee-9dad-a7934ea5d560
# ╠═59a80010-21a8-11f0-2aaa-5528a16a7081
# ╟─66ffc650-93e6-4705-8212-7958f75bb1e1
# ╟─e7657c9d-1ff1-46c8-a19e-99329474cae8
