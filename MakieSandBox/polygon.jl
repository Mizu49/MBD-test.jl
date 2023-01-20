using Revise
using GLMakie
using GeometryBasics

l1 = 2
phi1 = 7/4 * pi
l2 = 2
phi2 = 5/4 * pi

f = Figure()
ax = Axis(f[1, 1])

b1_p1 = Point(0.0, 0.0)
b1_p2 = Point(l1 * cos(phi1), l1 * sin(phi1));
b2_p1 = Point(l1 * cos(phi1), l1 * sin(phi1))
b2_p2 = Point(l1 * cos(phi1) + l2 * cos(phi2), l1 * sin(phi1) + l2 * sin(phi2));

poly!(Polygon([b1_p1, b1_p2]), color = :red, strokecolor = :red, strokewidth = 1)
poly!(Polygon([b2_p1, b2_p2]), color = :red, strokecolor = :blue, strokewidth = 1)

xlims!(-4, 4)
ylims!(-4, 0.5)
vlines!(ax, 0, color = :black)
hlines!(ax, 0, color = :black)
# hidespines!(ax)

display(f)
