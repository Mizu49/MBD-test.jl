using GeometryBasics, GLMakie, FileIO

# overload the multiplication operator 
Base.:*(C::AbstractMatrix, points::AbstractVector{<:AbstractPoint}) = map(point -> Point{3}(C * point) , points)

bunny = load("Stanford_Bunny.stl")

bunny_points = decompose(Point3f, bunny)
bunny_faces  = decompose(GLTriangleFace, bunny)

theta = pi/2
C = [
    1          0           0
    0 cos(theta) -sin(theta)
    0 sin(theta)  cos(theta)
]

# rotate all points with rotation matrix
bunny_points2 = C * bunny_points

fig = Figure(resolution = (600, 600))
ax = Axis3(
    fig[1, 1],
    aspect = :data,
    viewmode = :fit
)
mesh!(ax, bunny_points2, bunny_faces, color=:blue, linewidth=2)

display(fig)
