using Plots, LinearAlgebra

# パラメータ設定
l1 = 1.0
s1 = 0.5
g = 9.81
m1 = 5.0
I1 = m1 * l1^2 / 12

function EOM(time, state)
    
    x1      = state[1]
    y1      = state[2]
    phi1    = state[3]
    x1dot   = state[4]
    y1dot   = state[5]
    phi1dot = state[6]

    # 一般化質量行列
    M = diagm([m1, m1, I1])
    # 一般化外力ベクトル
    Q = [0, -m1*g, 0]
    # 拘束条件式
    C = [
        x1 - s1 * cos(phi1)
        y1 - s1 * sin(phi1)
    ]

    # ヤコビアン
    Cq = [
        1 0  s1 * sin(phi1)
        0 1 -s1 * cos(phi1)
    ]

    A = [
        M  transpose(Cq)
        Cq zeros(2, 2)
    ]
    
    Cdot = Cq * state[4:6]

    Gm = [
        -s1 * phi1dot^2 * cos(phi1)
        -s1 * phi1dot^2 * sin(phi1)
    ]

    # バウムガルテの安定化法
    alpha = 10
    beta = 10
    Gm = Gm - 2 * alpha * Cdot - beta^2 * C

    RHS = vcat(Q, Gm)

    accel_lambda  = A \ RHS

    # qdot qddot
    differential = vcat(state[4:6], accel_lambda[1:3])

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

timelength = 5.0
Ts = 1e-2

datanum = Integer(timelength / Ts + 1)
times = 0.0:Ts:timelength
states = [zeros(6) for _ in 1:datanum]
states[1] = [l1/2, 0.0, 0.0, 0.0, 0.0, 0.0]

for idx = 1:datanum-1
    
    states[idx+1] = runge_kutta(times[idx], states[idx], Ts)

end

plt = plot()
plot!(plt, times, getindex.(states, 3), label = "Angle")
xlabel!(plt, "Time (s)")
ylabel!(plt, "Angle (rad)")
display(plt)
