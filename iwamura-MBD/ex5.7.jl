using Plots, LinearAlgebra, Symbolics

@variables x1 y1 phi1 t

sym_q = [x1, y1, phi1]

C = [
    x1 - 0.8 * cos(phi1),
    y1 + 1.1 * sin(phi1),
    phi1 - 1/2 * t^2
]

Cq = Symbolics.jacobian(C, sym_q)
Ct = Symbolics.jacobian(C, [t])

funcC  = eval(build_function(C,  [t, x1, y1, phi1])[1])
funcCq = eval(build_function(Cq, [t, x1, y1, phi1])[1])
funcCt = eval(build_function(Ct, [t, x1, y1, phi1])[1])

initialval = [0.8, 0.0, 0.0]
timelength = 5.0
dt = 1e-2

datanum = Integer(timelength / dt + 1)

time = 0:dt:timelength
q_array = [zeros(3) for _ in 1:datanum]
q_array[1] = initialval
dq_array = [zeros(3) for _ in 1:datanum]

for idx = 1:datanum-1
    
    t = time[idx]

    # 位置解析
    iter = 0
    conv = 1
    while (iter < 100) && (conv == 1)
        
        augstate = vcat(t, q_array[idx])
        currentC = funcC(augstate)
        currentCq = funcCq(augstate)

        dq = - inv(currentCq) * currentC
        q_array[idx+1] = q_array[idx] + dq

        conv = 0
        for j = 1:3
            if (abs(currentC[j, 1]) > 0.0001)
                conv = 1
            end
        end
        iter = iter + 1
    end

    # 速度解析
    augstate = vcat(t, q_array[idx])
    currentCt = funcCt(augstate)
    currentCq = funcCq(augstate)

    dq_array[idx] = vec(inv(currentCq) * currentCt)

end

plotlyjs()
fig1 = plot()
plot!(fig1, time, getindex.(q_array, 1), label = "x1")
plot!(fig1, time, getindex.(q_array, 2), label = "y1")
plot!(fig1, time, getindex.(q_array, 3), label = "phi1")
