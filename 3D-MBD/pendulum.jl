using Revise, GLMakie, LinearAlgebra, StaticArrays, BlockDiagonals, GeometryBasics, Symbolics

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

@variables t, sym_q[1:7]
sym_q = collect(sym_q)

# オイラーパラメータ -> BRF
A_1 = [
    (1 -2 * sym_q[6]^2 -2 * sym_q[7]^2) (2 * (sym_q[5] * sym_q[6] - sym_q[4] * sym_q[7])) (2 * (sym_q[5] * sym_q[7] + sym_q[4] * sym_q[6]))
    (2 * (sym_q[5] * sym_q[6] + sym_q[4] * sym_q[7])) (1 -2 * sym_q[5]^2 -2 * sym_q[7]^2) (2 * (sym_q[6] * sym_q[7] - sym_q[4] * sym_q[5])) 
    (2 * (sym_q[4] * sym_q[6] - sym_q[4] * sym_q[6])) (2 * (sym_q[6] * sym_q[7] + sym_q[4] * sym_q[5])) (1 -2 * sym_q[4]^2 -2 * sym_q[5]^2)
]

# BRF -> point
u_bar_1 = [0, 0, s1]

# 拘束条件のシンボリック表現
C = [
    sym_q[1:3] - A_1 * u_bar_1
    transpose(sym_q[1:3]) * sym_q[1:3]  - s1^2
]

# ヤコビアンのシンボリック表現
Cq = Symbolics.jacobian(C, sym_q)
Ct = Symbolics.jacobian(C, [t])

# シンボリック表現を関数化
func_constraint = eval(build_function(C, sym_q)[1])
func_jacobian = eval(build_function(Cq, sym_q)[1])

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
    accel2 = [SVector{6}(zeros(6)) for _ in 1:datanum]

    print("Running simulation ...")
    simtime = @elapsed for idx = 1:datanum-1

        # DAEから加速度を計算
        accel[idx] = calc_acceleration(times[idx], states[idx])

        # 後退差分で加速度を計算
        if idx != 1
            accel2[idx] = (states[idx][7:12] - states[idx-1][7:12]) ./ Ts
        end

        # time evolution
        states[idx+1] = runge_kutta(times[idx], states[idx], Ts)

    end
    println("done!")

    fig1 = Figure()
    ax1 = Axis(fig1[1, 1])
    lines!(ax1, times, getindex.(states, 3))
    lines!(ax1, times, getindex.(states, 6))
    ax1.xlabel = "Time (s)"
    ax1.ylabel = "Angle (rad)"

    fig2 = Figure()
    ax2 = Axis(fig2[1, 1])
    lines!(ax2, times, getindex.(accel, 3))
    lines!(ax2, times, getindex.(accel, 6))
    lines!(ax2, times, getindex.(accel2, 3))
    lines!(ax2, times, getindex.(accel2, 6))
    ax2.xlabel = "Time (s)"
    ax2.ylabel = "Angular acceleration (rad/s^2)"

    figures = [fig1, fig2]

    print("Generating animation ...")
    anim_time = @elapsed generate_animation(times, states)
    println("done!")

    println("simulation: $simtime (s) / animation: $anim_time (s)")

    return (states, figures)
end

# (states, figures) = main();
