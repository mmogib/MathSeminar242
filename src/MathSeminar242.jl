### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 59a80010-21a8-11f0-2aaa-5528a16a7081
begin
	
	using Pkg: Pkg
	Pkg.activate("../")
	using Plots, PlotThemes
	using OAFCVX
	using CommonMark
	using PlutoUI
	using LaTeXStrings, Latexify, HypertextLiteral, Colors
	using LinearAlgebra
	md"``\color{white}{\text{Packages}}``"
end

# ╔═╡ 1c7f4f9c-0954-408f-927f-66d5929fc006
TableOfContents(title="📚 MATH SEMINAR: TERM 242", indent=true, depth=4)

# ╔═╡ 85b8d34e-d89d-49f0-b44e-fc0acee01d73
md"# Approximating Convex Vector Optimization Solutions"

# ╔═╡ 8c4c7c31-e590-4ff5-bc7b-f14d15293980
cm"""
This seminar provides an accessible introduction to scalarization techniques for approximating solutions of convex vector optimization problems. After motivating the general convex vector‐optimization framework and its “upper image,” we review key definitions—C-convexity, minimal and weakly minimal elements, polyhedral cones, and ε-solutions—through intuitive examples. 

Emphasis is placed on clear geometric interpretation and practical implementation, making these techniques readily approachable for newcomers. 

We conclude by sketching promising research directions: handling unbounded vector problems and establishing precise rates of convergence for these algorithms, with an eye toward both theoretical guarantees and computational efficiency.
"""

# ╔═╡ dc6dbd80-444f-4697-ad54-476b16cd5fe6
md"## The Problem"

# ╔═╡ 9f14a6ae-8756-4059-a984-3b935ff598a0
let
	# Pick a theme
	theme(:wong2)
	
	# Disk parameters
	radius = 1.0
	xc, yc = 1.0, 1.0
	
	# Parameterize the boundary circle
	θ = LinRange(0, 2π, 200)
	x = xc .+ radius .* cos.(θ)
	y = yc .+ radius .* sin.(θ)
	
	# Draw the filled disk
	plot(
	    x, y;
	    seriestype   = :shape,      # draw as a filled polygon
	    fillcolor    = :skyblue,    # fill color
	    fillalpha    = 0.4,         # transparency
	    linecolor    = :blue,       # boundary color
	    linewidth    = 0.5,
	    aspect_ratio = :equal,      # equal scaling on x and y
	    legend       = false,
	    xlabel       = "x",
	    ylabel       = "y",
	    title        = L"f(\Omega)=\Omega",
		frame_style  = :origin
	)
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
	# # Optionally mark the center
	# scatter!([xc], [yc];
	#     color   = :red,
	#     marker  = :circle,
	#     markersize = 6,
	#     label   = "Center"
	# )
end	

# ╔═╡ 764259e7-1d76-495b-8553-89f72caf3e46
cm"## Definitions and Background"

# ╔═╡ d9c32ddb-9922-4b73-882d-bc353e0dd761
let
	# pick a theme
	theme(:wong2)
	
	# 1) rectangle = positive‐orthant truncated to [0,2]×[0,2]
	rect_x1 = [1.0, 3.0, 3.0, 1.0]
	rect_y1 = [0.0, 0.0, 3.0, 3.0]
	rect_x2 = [0.0, 3.0, 3.0, 0.0]
	rect_y2 = [1.0, 1.0, 3.0, 3.0]
	
	# 2) circle = unit‐disk centered at (1,1)
	xc, yc, r = 1.0, 1.0, 1.0
	θ = range(0, 2π, length=300)
	circ_x = xc .+ r * cos.(θ)
	circ_y = yc .+ r * sin.(θ)
	
	# 3) start the plot by drawing the rectangle
	p = plot(
	    [(rect_x1, rect_y1),(rect_x2, rect_y2)],
	    seriestype   = :shape,
	    fillcolor    = :darkgray,
	    fillalpha    = 1,
	    linecolor    = :transparent,
	    label        = "",
	    xlims        = (-0.2, 3.2),
	    ylims        = (-0.2, 3.2),
	    aspect_ratio = :equal,
	    xlabel       = L"f_1",
	    ylabel       = L"f_2",
	    legend       = false,
	    grid         = false,
		frame_style  = :origin
	)
	
	# 4) overlay the filled disk
	plot!(
	    p, circ_x, circ_y,
	    seriestype   = :shape,
	    fillcolor    = :gray,
	    fillalpha    = 1,
	    linecolor    = :gray,
	    linewidth    = 1.5,
	    label        = "",
	)
	
	# 5) add the two text labels
	annotate!(p, 1.0, 1.0, text(L"f(\mathcal{\Omega})", :center, 14))
	annotate!(p, 1.8, 1.8, text(L"\mathcal{P}=f(\Omega)+\mathbb{R}^2_+",    :left,   12))	
end

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

# ╔═╡ a5bde7fe-472a-46e5-9cab-494e9c945493
let
	# pick a theme
	theme(:wong2)
	
	# 1) rectangle = positive‐orthant truncated to [0,2]×[0,2]
	rect_x1 = [1.0, 3.0, 3.0, 1.0]
	rect_y1 = [0.0, 0.0, 3.0, 3.0]
	rect_x2 = [0.0, 3.0, 3.0, 0.0]
	rect_y2 = [1.0, 1.0, 3.0, 3.0]
	
	# 2) circle = unit‐disk centered at (1,1)
	xc, yc, r = 1.0, 1.0, 1.0
	θ = range(0, 2π, length=300)
	circ_x = xc .+ r * cos.(θ)
	circ_y = yc .+ r * sin.(θ)
	
	# 3) start the plot by drawing the rectangle
	p = plot(
	    [(rect_x1, rect_y1),(rect_x2, rect_y2)],
	    seriestype   = :shape,
	    fillcolor    = :darkgray,
	    fillalpha    = 1,
	    linecolor    = :transparent,
	    label        = "",
	    xlims        = (-0.2, 3.2),
	    ylims        = (-0.2, 3.2),
	    aspect_ratio = :equal,
	    xlabel       = L"f_1",
	    ylabel       = L"f_2",
	    legend       = false,
	    grid         = false,
		frame_style  = :origin
	)
	
	# 4) overlay the filled disk
	plot!(
	    p, circ_x, circ_y,
	    seriestype   = :shape,
	    fillcolor    = :gray,
	    fillalpha    = 1,
	    linecolor    = :gray,
	    linewidth    = 1.5,
	    label        = "",
	)
	θ_t = π/4
	d = [cos(θ_t), sin(θ_t)]  # points toward the center
	quiver!(p, [0.0], [0.0],
	    quiver      = ([d[1]], [d[2]]),
	    arrow       = true,
	    linewidth   = 1.2,
	    linecolor   = :blue,
	    linestyle   = :dashdot,
	)

	scatter!(p, [0.0], [0.0];
	    color      = :black,
	    markersize = 4,
	)
	p_x = xc + r*cos(θ_t+π)
	p_y = yc + r*sin(θ_t+π)
	b = dot(d, [p_x, p_y])
	xs = [-0.2, 2.6]
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
	
	# 5) add the two text labels
	p_x = xc/8 + r*cos(θ_t)
	p_y = yc/12 + r*sin(θ_t)
	annotate!(p, 1.0, 1.2, text(L"f(\mathcal{\Omega})", :center, 14))
	annotate!(p, 2.8, 2.8, text(L"\mathcal{P}",    :left,   12))
	annotate!(p, p_x,  p_y,  text(L"d",             :blue, 14))
	annotate!(p, -0.3,  -.1,  text(L"v=\mathbf{0}",             14))
	annotate!(p, -0.3,  .61,  text(L"H",      :red,       14))
end

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

# ╔═╡ a707080d-33ef-4dda-b2a3-72a23dc90db3
let
	# pick a theme
	theme(:wong2)
	
	# 1) rectangle = positive‐orthant truncated to [0,2]×[0,2]
	rect_x1 = [1.0, 3.5, 3.5, 1.0]
	rect_y1 = [0.0, 0.0, 3.5, 3.5]
	rect_x2 = [0.0, 3.5, 3.5, 0.0]
	rect_y2 = [1.0, 1.0, 3.5, 3.5]
	
	# 2) circle = unit‐disk centered at (1,1)
	xc, yc, r = 1.0, 1.0, 1.0
	θ = range(0, 2π, length=300)
	circ_x = xc .+ r * cos.(θ)
	circ_y = yc .+ r * sin.(θ)
	
	# 3) start the plot by drawing the rectangle
	p = plot(
	    [(rect_x1, rect_y1),(rect_x2, rect_y2)],
	    seriestype   = :shape,
	    fillcolor    = :darkgray,
	    fillalpha    = 1,
	    linecolor    = :transparent,
	    label        = "",
	    xlims        = (-0.5, 3.5),
	    ylims        = (-0.5, 3.5),
	    aspect_ratio = :equal,
	    xlabel       = L"f_1",
	    ylabel       = L"f_2",
	    legend       = false,
	    grid         = false,
		frame_style  = :origin
	)
	
	# 4) overlay the filled disk
	plot!(
	    p, circ_x, circ_y,
	    seriestype   = :shape,
	    fillcolor    = :gray,
	    fillalpha    = 1,
	    linecolor    = :gray,
	    linewidth    = 1.5,
	    label        = "",
	)

	# 5) add the two text labels
	annotate!(p, 1.0, 1.0, text(L"f(\mathcal{\Omega})", :center, 14))
	
	# 6) inner P
	xs = [10.0, 1.0, 1-1/sqrt(2), 0.0, 0.0]
	ys = [0.0, 0.0, 1-1/sqrt(2), 1.0, 10.0]
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
	annotate!(p, 1.2, 2.8, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+)",    :left,   12))	
	
	# 7) inner P
	xs = [10.0, 1.0, 1-1/sqrt(2), 0.0, 0.0]
	ys = [0.0, 0.0, 1-1/sqrt(2), 1.0, 10.0]
	ϵ = 0.06
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
	annotate!(p, 1.2, 2.8, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+)",    :left,   12,"#483D8B"))	
	annotate!(p, 1.2, 2.4, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+ -\mathbb{1}\{\epsilon\})",    :left,   12,dark_palette[:dark_magenta]))	
end


# ╔═╡ 5b0e766a-28fd-467d-8c70-717730ec47ef
let
	# pick a theme
	theme(:wong2)
	
	# 1) rectangle = positive‐orthant truncated to [0,2]×[0,2]
	rect_x1 = [1.0, 3.5, 3.5, 1.0]
	rect_y1 = [0.0, 0.0, 3.5, 3.5]
	rect_x2 = [0.0, 3.5, 3.5, 0.0]
	rect_y2 = [1.0, 1.0, 3.5, 3.5]
	
	# 2) circle = unit‐disk centered at (1,1)
	xc, yc, r = 1.0, 1.0, 1.0
	θ = range(0, 2π, length=300)
	circ_x = xc .+ r * cos.(θ)
	circ_y = yc .+ r * sin.(θ)
	
	# 3) start the plot by drawing the rectangle
	p = plot(
	    [(rect_x1, rect_y1),(rect_x2, rect_y2)],
	    seriestype   = :shape,
	    fillcolor    = :darkgray,
	    fillalpha    = 1,
	    linecolor    = :transparent,
	    label        = "",
	    xlims        = (-0.5, 3.5),
	    ylims        = (-0.5, 3.5),
	    aspect_ratio = :equal,
	    xlabel       = L"f_1",
	    ylabel       = L"f_2",
	    legend       = false,
	    grid         = false,
		frame_style  = :origin
	)
	
	# 4) overlay the filled disk
	plot!(
	    p, circ_x, circ_y,
	    seriestype   = :shape,
	    fillcolor    = :gray,
	    fillalpha    = 1,
	    linecolor    = :gray,
	    linewidth    = 1.5,
	    label        = "",
	)

	# 5) add the two text labels
	annotate!(p, 1.0, 1.0, text(L"f(\mathcal{\Omega})", :center, 14))
	
	# 6) inner P
	xs = [10.0, 1.0, 1-1/sqrt(2), 0.0, 0.0]
	ys = [0.0, 0.0, 1-1/sqrt(2), 1.0, 10.0]
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
	annotate!(p, 1.2, 2.8, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+)",    :left,   12))	
	
	# 7) inner P
	xs = [10.0, 1.0, 1-1/sqrt(2), 0.0, 0.0]
	ys = [0.0, 0.0, 1-1/sqrt(2), 1.0, 10.0]
	ϵ = 0.06
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
	annotate!(p, 1.2, 2.8, text(L"\operatorname{bd}(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+)",    :left,   12,"#483D8B"))	
	annotate!(p, 1.2, 2.4, text(L"\operatorname{bd}\left(\operatorname{conv} f(\overline{\Omega})+\mathbb{R}^2_+ +\mathbf{B}(\mathbf{0},\epsilon)\right)",    :left,   12,dark_palette[:dark_magenta]))	
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
&  \|x-\mathbb{1}\|_2^2 \leq 1
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
C^{+}:=\operatorname{coconv}\left\{w^1, \ldots, w^l\right\}, 
``` 
and ``\mathcal{\Omega}`` is compact with nonempty interior.
$(ebl())

$(define("Polyhedral Convex Set"))
A set ``S\in\mathbb{R}^m`` is a __polyhedral convex set__ if it is of the form 
```math
S=\left\{y \in \mathbb{R}^m \mid A^T y \geq b\right\},
```
where ``A \in`` ``\mathbb{R}^{m \times n}, b \in \mathbb{R}^n``. 
$(ebl())

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
\mathbb{R}^m_+=\left\{y \in \mathbb{R}^m \mid I y \geq \mathbf{0}\right\},
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
# ╟─57b118d7-a7dd-4f57-95e3-ef7d3fc334bf
# ╟─9f14a6ae-8756-4059-a984-3b935ff598a0
# ╟─764259e7-1d76-495b-8553-89f72caf3e46
# ╟─8a7307cb-0ad4-4e46-b18d-92eb89ac82cb
# ╟─871cf958-6093-4004-99fd-ead3a2561e50
# ╟─647e1f03-9a12-4239-b649-8d9aabbaa7c2
# ╟─3bf3f0c2-496f-41b5-b2bc-0c0cb7b6feab
# ╟─d9c32ddb-9922-4b73-882d-bc353e0dd761
# ╟─bb397f1a-1cde-4e83-afe3-e9dbad4a485f
# ╟─035fe90d-6c65-4a6a-b345-10d44e12b131
# ╟─a707080d-33ef-4dda-b2a3-72a23dc90db3
# ╟─58d3a13b-87cd-481c-a208-7259607b034c
# ╟─5b0e766a-28fd-467d-8c70-717730ec47ef
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
# ╟─a5bde7fe-472a-46e5-9cab-494e9c945493
# ╟─f9a1dfe1-bc71-42c4-a6e8-3f422eba7507
# ╟─e3f5835c-18fc-4f00-92f5-e55598e68081
# ╠═c6db5f72-a8fa-481c-85e9-7c587e0589b4
# ╠═34c3882f-2178-40f8-80b8-e3f019e5284e
# ╟─59a80010-21a8-11f0-2aaa-5528a16a7081
# ╟─66ffc650-93e6-4705-8212-7958f75bb1e1
# ╟─e7657c9d-1ff1-46c8-a19e-99329474cae8
