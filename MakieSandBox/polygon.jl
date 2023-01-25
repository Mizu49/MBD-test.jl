using GeometryBasics, GLMakie, FileIO

bunny = load("Stanford_Bunny.stl")

@time begin
    f1 = Figure(resolution = (600, 600))
    ax1 = Axis3(f1[1, 1]) # 3D plotする時はAxis3で3次元の軸を作る
    poly!(ax1, bunny, color=:blue, linewidth=2)    
end

@time begin
    f2 = Figure(resolution = (600, 600))
    ax2 = Axis3(f2[1, 1]) # 3D plotする時はAxis3で3次元の軸を作る
    mesh!(ax2, bunny, color=:blue, linewidth=2)
end

display(f1)
display(f2)
