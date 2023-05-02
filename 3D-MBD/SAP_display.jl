using GLMakie, GeometryBasics

struct PolygonPlate
    a::Real # x-direction length
    b::Real # y-direction length
    t::Real # z-direction length (thickness)
end


function polygon_plate(pos::AbstractVector, p::PolygonPlate)
   
    C = quaternion2dcm(pos[4:7])

    points = hcat(
        pos[1:3] + C * [-p.a/2, -p.b/2, +p.t/2],
        pos[1:3] + C * [-p.a/2, +p.b/2, +p.t/2],
        pos[1:3] + C * [+p.a/2, +p.b/2, +p.t/2],
        pos[1:3] + C * [+p.a/2, -p.b/2, +p.t/2],
        pos[1:3] + C * [-p.a/2, -p.b/2, -p.t/2],
        pos[1:3] + C * [-p.a/2, +p.b/2, -p.t/2],
        pos[1:3] + C * [+p.a/2, +p.b/2, -p.t/2],
        pos[1:3] + C * [+p.a/2, -p.b/2, -p.t/2]
    )
    
    faces = [
        1 2 3
        3 4 1
        1 2 6
        1 5 6
        1 4 5
        4 5 8
        3 4 7
        4 7 8
        5 6 7
        5 7 8
        2 3 6
        3 6 7
    ]
    
    return (points, faces)
end

params = PolygonPlate(l1, l1, 0.01)

q = states[5000]
(points, faces) = polygon_plate(q, params)

obs_points = Observable(points)

fig = Figure()
ax = Axis3(
    fig[1, 1],
    # xlabel = "x label",
    # ylabel = "y label",
    # zlabel = "z label",
    elevation = pi/6,
    azimuth   = pi/4,
    aspect = :data,
    viewmode = :fit,
    limits = (-1.2, 1.2, -1.2, 1.2, -1.2, 1.2)
)
hidespines!(ax)
hidedecorations!(ax)

# Axis
lines!(ax, [-1.2, 1.2], [0.0, 0.0], [0.0, 0.0], linestyle = :dash, color = :black, linewidth = 1)
lines!(ax, [0.0, 0.0], [-1.2, 1.2], [0.0, 0.0], linestyle = :dash, color = :black, linewidth = 1)
lines!(ax, [0.0, 0.0], [0.0, 0.0], [-1.2, 1.2], linestyle = :dash, color = :black, linewidth = 1)

# 平板のプロットを実行
mesh!(ax, obs_points, faces, color=:blue, shading = true)

iter = 1:100:size(states, 1)
record(fig, "animation.gif", iter, framerate = 15) do idx

    q = states[idx]
    (points, faces) = polygon_plate(q, params)

    obs_points[] = points
end
