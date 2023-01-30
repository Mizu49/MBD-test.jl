using GeometryBasics, GLMakie, FileIO

# overload the multiplication operator 
Base.:*(C::AbstractMatrix, points::AbstractVector{<:AbstractPoint{Dim, T}}) where {Dim, T} = map(point -> Point{Dim, T}(C * point) , points)

bunny = load("Stanford_Bunny.stl")

bunny_points = decompose(Point3f, bunny)
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

fig = Figure(resolution = (600, 600))
ax = Axis3(
    fig[1, 1],
    aspect = :data,
    viewmode = :fit,
    limits = (-100, 100, -100, 100, 0, 100)
)
mesh!(ax, points, bunny_faces, color=:blue, linewidth=2)

ratio = 2 * pi / 60

iter = 0:120
record(fig, "animation.gif", iter; framerate = 20) do idx
    points[] = rotate_bunny(points[], ratio)
end