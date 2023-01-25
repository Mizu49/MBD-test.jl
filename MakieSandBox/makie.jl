using Revise, GLMakie, GeometryBasics

# true: show in VScode plot navigator
Makie.inline!(false)

obj = Rect3(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)

fig1 = Figure(resolution = (600, 600))
ax1 = Axis3( # 3D plotする時はAxis3で3次元の軸を作る
    fig1[1, 1],
    title = "mesh plot",
    xlabel = "x label",
    ylabel = "y label",
    zlabel = "z label",
    elevation = pi/4,
    azimuth   = pi/4,
    aspect = :data,
    viewmode = :fit,
) 

hidespines!(ax1)
hidedecorations!(ax1)

mesh!(ax1, obj, color=:blue, linewidth = 2)

display(fig1)
