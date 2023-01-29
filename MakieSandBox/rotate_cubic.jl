using GLMakie, GeometryBasics

theta = pi/4
C = [
    1          0           0
    0 cos(theta) -sin(theta)
    0 sin(theta)  cos(theta)
]

points = [
    Point3{Float64}(C * [ 0.5,  0.5,  0.5]),
    Point3{Float64}(C * [-0.5,  0.5,  0.5]),
    Point3{Float64}(C * [-0.5, -0.5,  0.5]),
    Point3{Float64}(C * [ 0.5, -0.5,  0.5]),
    Point3{Float64}(C * [ 0.5,  0.5, -0.5]),
    Point3{Float64}(C * [-0.5,  0.5, -0.5]),
    Point3{Float64}(C * [-0.5, -0.5, -0.5]),
    Point3{Float64}(C * [ 0.5, -0.5, -0.5])
]

cube = GeometryBasics.mesh(points)

fig = Figure()
ax = Axis3(
    fig[1, 1],
    # elevation = pi/6,
    # azimuth   = pi/4,
    aspect = :data,
    viewmode = :fit
    ) # 3D plotする時はAxis3で3次元の軸を作る
mesh!(ax, cube, color=:yellow, shading = true)

display(fig)

