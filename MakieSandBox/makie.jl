using GLMakie, Plots

fig1 = Plots.plot()
gr()
Plots.plot!(fig1, rand(10), rand(10))

display(fig1)

fig2 = Figure()
Axis(fig2[1, 1])
GLMakie.plot!(rand(10), rand(10))
