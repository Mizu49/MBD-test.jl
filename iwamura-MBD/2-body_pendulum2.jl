using Revise, Plots, LinearAlgebra, StaticArrays, BlockDiagonals

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

# エレメントの質量行列
function M_elem(mass, inertia)
    M = SMatrix{3, 3}(diagm([mass, mass, inertia]))
    return M
end

# システムの質量行列
function func_global_mass(element_masses::Vector{<:AbstractMatrix})
    # 要素数
    elemnum = size(element_masses, 1)
    
    # ブロック対角行列にする
    return SMatrix{3 * elemnum, 3 * elemnum}(BlockDiagonal(element_masses))
end

M1 = M_elem(m1, I1)
M2 = M_elem(m2, I2)

"""
拘束条件式
"""
function func_constraint(q::SVector{6})::SVector
    
    C = SVector{4}([
        q[1] - s1 * cos(q[3])
        q[2] - s1 * sin(q[3])
        q[1] + (l1 - s1) * cos(q[3]) - q[4] + s2 * cos(q[6])
        q[2] + (l1 - s1) * sin(q[3]) - q[5] + s2 * sin(q[6])
    ])

    return C
end

function func_jacobian(q::SVector{6})::SMatrix
    
    # ヤコビアン
    Cq = SMatrix{4, 6}([
        1 0  s1 * sin(q[3])        0   0   0
        0 1 -s1 * cos(q[3])        0   0   0
        1 0 -(l1 - s1) * sin(q[3]) -1  0   -s2 * sin(q[6])
        0 1 (l1 - s1) * cos(q[3])  0   -1  s2 * cos(q[6])
    ])

    return Cq
end

"""
ベクトル γ
"""
function func_gamma(q, qdot)
    
    gamma = [
        -s1 * qdot[3]^2 * cos(q[3])
        -s1 * qdot[3]^2 * sin(q[3])
        (l1 - s1) * qdot[3]^2 * cos(q[3]) + s2 * qdot[6]^2 * cos(q[6])
        (l1 - s1) * qdot[3]^2 * sin(q[3]) + s2 * qdot[6]^2 * sin(q[6])
    ]

    return gamma
end

"""
一般化外力
"""
function func_external_force()::SVector

    Q = SVector{6}([0, -m1*g, 0, 0, -m2*g, 0])
    
    return Q
end

function EOM(time, state)
    
    # 状態量をパースする
    q    = SVector{6}(state[1:6])
    qdot = SVector{6}(state[7:12])

    # 一般化質量行列
    M = func_global_mass([M1, M2])

    # 一般化外力ベクトル
    Q = func_external_force()

    # 拘束条件式
    C = func_constraint(q)
    
    # ヤコビアン
    Cq = func_jacobian(q)
    
    Cdot = Cq * qdot

    Gm = func_gamma(q, qdot)

    A = [
        M  transpose(Cq)
        Cq zeros(4, 4)
    ]
    
    # バウムガルテの安定化法
    alpha = 10
    beta = 10
    Gm = Gm - 2 * alpha * Cdot - beta^2 * C

    RHS = vcat(Q, Gm)

    accel_lambda  = A \ RHS

    # qdot qddot
    differential = vcat(qdot, accel_lambda[1:6])

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

function main()
    timelength = 20.0
    Ts = 1e-2

    datanum = Integer(timelength / Ts + 1)
    times = 0.0:Ts:timelength
    states = [SVector{12}(zeros(6*2)) for _ in 1:datanum]
    states[1] = vcat([l1/2, 0.0, 0.0, l1 + l2/2, 0.0, 0.0], zeros(6))
    accel = [SVector{6}(zeros(6)) for _ in 1:datanum]

    for idx = 1:datanum-1

        gamma = func_gamma(states[idx][1:6], states[idx][7:12])
        Cq = func_jacobian(SVector{6}(states[idx][1:6]))

        # ここは疑似逆行列でいいのか？？？
        accel[idx] = pinv(Cq) * gamma

        # time evolution
        states[idx+1] = runge_kutta(times[idx], states[idx], Ts)

    end

    fig1 = plot()
    plot!(fig1, times, getindex.(states, 3), label = "phi1")
    plot!(fig1, times, getindex.(states, 6), label = "phi2")
    xlabel!(fig1, "Time (s)")
    ylabel!(fig1, "Angle (rad)")
    title!(fig1, "Angle of joint")
    display(fig1)

    fig2 = plot()
    plot!(fig2, times, getindex.(accel, 3), label = "phi1 accel")
    plot!(fig2, times, getindex.(accel, 6), label = "phi2 accel")
    xlabel!(fig2, "Time (s)")
    ylabel!(fig2, "Angular acceleration (rad/s^2)")
    title!(fig2, "Angular acceleration")
    display(fig2)
end

main()
