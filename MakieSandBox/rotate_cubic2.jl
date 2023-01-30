using GLMakie, GeometryBasics

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

function rotate_cube(theta)
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
end

points = Observable(rotate_cube(0))

fig = Figure()
ax = Axis3(
    fig[1, 1],
    # elevation = pi/6,
    # azimuth   = pi/4,
    aspect = :data,
    viewmode = :fit
    )
mesh!(ax, points, faces, color = :blue ,shading = true)

hidedecorations!(ax)
hidespines!(ax)

ratio = pi/60

iter = 0:60
record(fig, "animation.gif", iter; framerate = 30) do idx
    points[] = rotate_cube(ratio * idx)
end