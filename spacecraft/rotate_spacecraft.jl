using GLMakie, GeometryBasics, StaticArrays

include("rotation_utilities.jl")

include("spacecraft_body.jl")

hoge = rand(3, 3)

spacecraft = get_spacecraft_polygon()

points = Observable(spacecraft.points)

fig = Figure()
ax = Axis3(
    fig[1, 1],
    # elevation = pi/6,
    # azimuth   = pi/4,
    aspect = :data,
    viewmode = :fit
    )
mesh!(ax, points, spacecraft.faces, color = :yellow ,shading = true)

hidedecorations!(ax)
hidespines!(ax)

ratios = SVector{3}([pi/30, pi/60, 0])

iter = 0:60
record(fig, "spacecraft.gif", iter; framerate = 30) do idx
    angles = ratios .* idx
    C = euler2dcm(angles)
    points[] = C * spacecraft.points
end