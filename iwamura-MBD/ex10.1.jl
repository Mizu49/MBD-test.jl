using Revise, Plots, LinearAlgebra, StaticArrays

# パラメータ設定
l1 = 2.0
s1 = l1/2
m1 = 0.5
I1 = m1 * l1^2 / 12

l2 = 2.0
s2 = l2/2
m2 = 0.5
I2 = m2 * l2^2 / 12

g = 9.81

function EOM(time, state)
    
    x1      = state[1]
    y1      = state[2]
    phi1    = state[3]
    x2      = state[4]
    y2      = state[5]
    phi2    = state[6]
    x1dot   = state[7]
    y1dot   = state[8]
    phi1dot = state[9]
    x2dot   = state[10]
    y2dot   = state[11]
    phi2dot = state[12]

    # 一般化質量行列
    M = diagm([m1, m1, I1, m2, m2, I2])
    # 一般化外力ベクトル
    Q = [0, -m1*g, 0, 0, -m2*g, 0]
    # 拘束条件式
    C = [
        x1 - s1 * cos(phi1)
        y1 - s1 * sin(phi1)
        x1 + (l1 - s1) * cos(phi1) - x2 + s2 * cos(phi2)
        y1 + (l1 - s1) * sin(phi1) - y2 + s2 * sin(phi2)
    ]

    # ヤコビアン
    Cq = [
        1 0  s1 * sin(phi1)        0   0   0
        0 1 -s1 * cos(phi1)        0   0   0
        1 0 -(l1 - s1) * sin(phi1) -1  0   -s2 * sin(phi2)
        0 1 (l1 - s1) * cos(phi1)  0   -1  s2 * cos(phi2)
    ]

    A = [
        M  transpose(Cq)
        Cq zeros(4, 4)
    ]
    
    Cdot = Cq * state[7:12]

    Gm = [
        -s1 * phi1dot^2 * cos(phi1)
        -s1 * phi1dot^2 * sin(phi1)
        (l1 - s1) * phi1dot^2 * cos(phi1) + s2 * phi2dot^2 * cos(phi2)
        (l1 - s1) * phi1dot^2 * sin(phi1) + s2 * phi2dot^2 * sin(phi2)
    ]

    # バウムガルテの安定化法
    alpha = 10
    beta = 10
    Gm = Gm - 2 * alpha * Cdot - beta^2 * C

    RHS = vcat(Q, Gm)

    accel_lambda  = A \ RHS

    # qdot qddot
    differential = vcat(state[7:12], accel_lambda[1:6])

    return differential
end

function runge_kutta(time, state, Ts)
    
    k1 = EOM(time, state)
    k2 = EOM(time + Ts/2, state + Ts/2 * k1)
    k3 = EOM(time + Ts/2, state + Ts/2 * k2)
    k4 = EOM(time + Ts  , state + Ts   * k3)

    nextstate = state + Ts/6 * (k1 + 2 * k2 + 2 * k3 + k4)

    return nextstate
end

timelength = 20.0
Ts = 1e-2

datanum = Integer(timelength / Ts + 1)
times = 0.0:Ts:timelength
states = [SVector{12}(zeros(6*2)) for _ in 1:datanum]
states[1] = vcat([l1/2, 0.0, 0.0, l1 + l2/2, 0.0, 0.0], zeros(6))
accel = [SVector{6}(zeros(6)) for _ in 1:datanum]

for idx = 1:datanum-1
    
    accel[idx] = states[idx][7:12]

    # time evolution
    states[idx+1] = runge_kutta(times[idx], states[idx], Ts)

end

plt = plot()
plot!(plt, times, getindex.(states, 3), label = "phi1")
plot!(plt, times, getindex.(states, 6), label = "phi2")
xlabel!(plt, "Time (s)")
ylabel!(plt, "Angle (rad)")
display(plt)
