using GeometryBasics, GLMakie


polygon = Rect3(0.0, 0.0, 0.0, 0.75, 0.5, 0.5)

f1 = Figure()
ax1 = Axis3(f1[1, 1]) # 3D plotする時はAxis3で3次元の軸を作る
poly!(ax1, polygon, color=:blue, linewidth=2)

f2 = Figure()
ax2 = Axis3(f2[1, 1]) # 3D plotする時はAxis3で3次元の軸を作る
mesh!(ax2, polygon, color=:blue, linewidth=2)

display(f1)
display(f2)
