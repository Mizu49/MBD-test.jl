using GLMakie, GeometryBasics, StaticArrays

include("rotation_utilities.jl")

include("spacecraft_body.jl")

label_positions = [
    Point3{Float64}(1.0, 0.0, 0.0),
    Point3{Float64}(-1.0, 0.0, 0.0),
    Point3{Float64}(0.0, 1.0, 0.0),
    Point3{Float64}(0.0, -1.0, 0.0),
    Point3{Float64}(0.0, 0.0, 1.0),
    Point3{Float64}(0.0, 0.0, -1.0),
]

label_texts = [
    "FWD",
    "AFT",
    "PRT",
    "STB",
    "ZNT",
    "NDR"
]

spacecraft = get_spacecraft_polygon()

points = Observable(spacecraft.points)

fig = Figure()
ax = Axis3(
    fig[1, 1],
    elevation = pi/6,
    azimuth   = pi/6,
    aspect = :data,
    viewmode = :fit,
    limits = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
)
mesh!(ax, points, spacecraft.faces, color = :yellow ,shading = true)

text!(
    label_positions,
    text = label_texts,
    align = (:center, :center),
)

lines!(ax, [-1, 1], [0, 0], [0, 0], linestyle = :dash, color = :black, linewidth = 3)
lines!(ax, [0, 0], [-1, 1], [0, 0], linestyle = :dash, color = :black, linewidth = 3)
lines!(ax, [0, 0], [0, 0], [-1, 1], linestyle = :dash, color = :black, linewidth = 3)

# hidedecorations!(ax)
# hidespines!(ax)

ratios = SVector{3}([pi/30, 0, 0])

iter = 0:300
record(fig, "spacecraft.gif", iter; framerate = 30) do idx    
    angles = ratios .* idx
    C = euler2dcm(angles)
    points[] = C * spacecraft.points
end