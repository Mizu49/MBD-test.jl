using GLMakie, GeometryBasics

struct PolygonPlate
    a::Real # x-direction length
    b::Real # y-direction length
    t::Real # z-direction length (thickness)
end


function polygon_plate(pos::AbstractVector, p::PolygonPlate)
   

    println(pos[1:3])

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

q = states[10000]

(points, faces) = polygon_plate(q, params)

fig = Figure()
ax = Axis3(
    fig[1, 1],
    # elevation = pi/6,
    # azimuth   = pi/4,
    aspect = :data,
    viewmode = :fit
    )
mesh!(ax, points, faces, color=:blue, shading = true)

