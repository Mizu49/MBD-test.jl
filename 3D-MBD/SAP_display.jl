using GLMakie, GeometryBasics, StaticArrays

struct PlateParameters
    a::Real # x-direction length
    b::Real # y-direction length
    t::Real # z-direction length (thickness)
end

struct PlatePolygon
    points::SizedMatrix # mutable statically-sized array
    faces::SMatrix # immutable statically-sized array
    parameters::PlateParameters
end

function PlatePolygon(coordinate::AbstractVector, parameters::PlateParameters)
    
    C = quaternion2dcm(coordinate[4:7])

    points = SizedMatrix{3, 8}(hcat(
        coordinate[1:3] + C * [-parameters.a/2, -parameters.b/2, +parameters.t/2],
        coordinate[1:3] + C * [-parameters.a/2, +parameters.b/2, +parameters.t/2],
        coordinate[1:3] + C * [+parameters.a/2, +parameters.b/2, +parameters.t/2],
        coordinate[1:3] + C * [+parameters.a/2, -parameters.b/2, +parameters.t/2],
        coordinate[1:3] + C * [-parameters.a/2, -parameters.b/2, -parameters.t/2],
        coordinate[1:3] + C * [-parameters.a/2, +parameters.b/2, -parameters.t/2],
        coordinate[1:3] + C * [+parameters.a/2, +parameters.b/2, -parameters.t/2],
        coordinate[1:3] + C * [+parameters.a/2, -parameters.b/2, -parameters.t/2]
    ))
    
    faces = SMatrix{12, 3}([
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
    ])
    
    return PlatePolygon(points, faces, parameters)
end

function update_polygon!(polygon::PlatePolygon, coordinate::AbstractVector)

    C = quaternion2dcm(coordinate[4:7])

    points = hcat(
        coordinate[1:3] + C * [-polygon.parameters.a/2, -polygon.parameters.b/2, +polygon.parameters.t/2],
        coordinate[1:3] + C * [-polygon.parameters.a/2, +polygon.parameters.b/2, +polygon.parameters.t/2],
        coordinate[1:3] + C * [+polygon.parameters.a/2, +polygon.parameters.b/2, +polygon.parameters.t/2],
        coordinate[1:3] + C * [+polygon.parameters.a/2, -polygon.parameters.b/2, +polygon.parameters.t/2],
        coordinate[1:3] + C * [-polygon.parameters.a/2, -polygon.parameters.b/2, -polygon.parameters.t/2],
        coordinate[1:3] + C * [-polygon.parameters.a/2, +polygon.parameters.b/2, -polygon.parameters.t/2],
        coordinate[1:3] + C * [+polygon.parameters.a/2, +polygon.parameters.b/2, -polygon.parameters.t/2],
        coordinate[1:3] + C * [+polygon.parameters.a/2, -polygon.parameters.b/2, -polygon.parameters.t/2]
    )

    # modify points of polygons
    polygon.points[:] = points

    return    
end

# specify the dimension of the plates
thickness = 0.05
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
    update_polygon!(plate1_polygon, q)

    ax.title = "Time: $(times[idx]) (s)"

    plate1_obs_points[] = plate1_polygon.points
end
