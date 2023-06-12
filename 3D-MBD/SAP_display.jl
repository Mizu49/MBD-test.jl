using GLMakie, GeometryBasics

struct PlateParameters
    a::Real # x-direction length
    b::Real # y-direction length
    t::Real # z-direction length (thickness)
end

struct PlatePolygon
    points
    faces
end


function PlatePolygon(coordinate::AbstractVector, p::PlateParameters)
   
    C = quaternion2dcm(coordinate[4:7])

    points = hcat(
        coordinate[1:3] + C * [-p.a/2, -p.b/2, +p.t/2],
        coordinate[1:3] + C * [-p.a/2, +p.b/2, +p.t/2],
        coordinate[1:3] + C * [+p.a/2, +p.b/2, +p.t/2],
        coordinate[1:3] + C * [+p.a/2, -p.b/2, +p.t/2],
        coordinate[1:3] + C * [-p.a/2, -p.b/2, -p.t/2],
        coordinate[1:3] + C * [-p.a/2, +p.b/2, -p.t/2],
        coordinate[1:3] + C * [+p.a/2, +p.b/2, -p.t/2],
        coordinate[1:3] + C * [+p.a/2, -p.b/2, -p.t/2]
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
    
    return PlatePolygon(points, faces)
end

function update_polygon!(polygon::PlatePolygon, p::PlateParameters, coordinate::AbstractVector)

    C = quaternion2dcm(coordinate[4:7])

    points = hcat(
        coordinate[1:3] + C * [-p.a/2, -p.b/2, +p.t/2],
        coordinate[1:3] + C * [-p.a/2, +p.b/2, +p.t/2],
        coordinate[1:3] + C * [+p.a/2, +p.b/2, +p.t/2],
        coordinate[1:3] + C * [+p.a/2, -p.b/2, +p.t/2],
        coordinate[1:3] + C * [-p.a/2, -p.b/2, -p.t/2],
        coordinate[1:3] + C * [-p.a/2, +p.b/2, -p.t/2],
        coordinate[1:3] + C * [+p.a/2, +p.b/2, -p.t/2],
        coordinate[1:3] + C * [+p.a/2, -p.b/2, -p.t/2]
    )

    # modify array
    polygon.points[:] = points

    return    
end

# specify the dimension of the plates
thickness = 0.01
plate1_params = PlateParameters(l1, l1, thickness)

# create inital polygon
q = states[1]
plate1_polygon = PlatePolygon(q, plate1_params)

# make observables
plate1_obs_points = Observable(plate1_polygon.points)

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
mesh!(ax, plate1_obs_points, plate1_polygon.faces, color=:blue, shading = true)

iter = 1:100:size(states, 1)
record(fig, "animation.gif", iter, framerate = 15) do idx

    q = states[idx]
    update_polygon!(plate1_polygon, plate1_params, q)

    ax.title = "Time: $(times[idx]) (s)"

    plate1_obs_points[] = plate1_polygon.points
end
