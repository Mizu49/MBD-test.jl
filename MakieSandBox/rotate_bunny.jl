using GeometryBasics, GLMakie, FileIO, Statistics

# overload the multiplication operator 
Base.:*(C::AbstractMatrix, points::AbstractVector{<:AbstractPoint{Dim, T}}) where {Dim, T} = map(point -> Point{Dim, T}(C * point) , points)

bunny = load("Stanford_Bunny.stl")

bunny_points = decompose(Point3f, bunny)
bunny_points = bunny_points .- mean(bunny_points)

bunny_faces  = decompose(GLTriangleFace, bunny)

function rotate_bunny(points, theta)

    C = [
        cos(theta) -sin(theta) 0
        sin(theta)  cos(theta) 0
                 0           0 1
    ]
    
    # rotate all points with rotation matrix
    return C * points
end

points = Observable(bunny_points)

slim = 60.0

xy_plane = Point{3, Float32}[(slim, slim, 0), (-slim, slim, 0), (-slim, -slim, 0), (slim, -slim, 0)]
yz_plane = Point{3, Float32}[(0, slim, slim), (0, -slim, slim), (0, -slim, -slim), (0, slim, -slim)]
zx_plane = Point{3, Float32}[(slim, 0, slim), (-slim, 0, slim), (-slim, 0, -slim), (slim, 0, -slim)]
plane = [
    1 2 3
    1 3 4
]

fig = Figure(resolution = (600, 600))
ax = Axis3(
    fig[1, 1],
    aspect = :data,
    viewmode = :fit,
    limits = (-slim, slim, -slim, slim, -slim, slim)
)
mesh!(ax, points, bunny_faces, color=:blue, linewidth=2)

poly!(ax, xy_plane, plane, color = (:gray, 0.2))
poly!(ax, yz_plane, plane, color = (:gray, 0.2))
poly!(ax, zx_plane, plane, color = (:gray, 0.2))

hidedecorations!(ax)
hidespines!(ax)

ratio = 2 * pi / 60

iter = 0:120
record(fig, "animation.gif", iter; framerate = 20) do idx
    points[] = rotate_bunny(points[], ratio)
end