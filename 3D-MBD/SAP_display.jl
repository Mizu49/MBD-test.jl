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
    
    C = transpose(quaternion2dcm(coordinate[4:7]))

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

    C = transpose(quaternion2dcm(coordinate[4:7]))

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

# numbers of the bodies
num_bodies = 2

# specify the dimension of the plates
thickness = 0.05
plate_params = [
    PlateParameters(l1, l1, thickness)
    PlateParameters(l2, l2, thickness)
]

# create inital polygon
currentstate = states[1]
plate_polygons = [
    PlatePolygon(currentstate[(1+7*idx):(7+7*idx)], plate_params[idx+1]) for idx in 0:num_bodies-1
]

# make observables
plate_obs_points = [
    Observable(plate_polygons[idx].points) for idx in 1:num_bodies
]

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
for idx = 1:num_bodies 
    mesh!(ax, plate_obs_points[idx], plate_polygons[idx].faces, color=(:blue, 0.3), shading = true)
end

iter = 1:100:size(states, 1)
record(fig, "animation.gif", iter, framerate = 15) do idx

    currentstate = states[idx]
    for idx = 0:num_bodies-1
        update_polygon!(plate_polygons[idx+1], currentstate[(1+7*idx):(7+7*idx)])
        plate_obs_points[idx+1][] = plate_polygons[idx+1].points
    end

    ax.title = "Time: $(times[idx]) (s)"

end
