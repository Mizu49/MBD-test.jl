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

dim_coordinate = 7

@variables t, sym_q[1:7], sym_qdot[1:7]
sym_q = collect(sym_q)
sym_qdot = collect(sym_qdot)

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

# 時間微分オペレータ
D = Differential(t)

# ヤコビアンのシンボリック表現
Cq = Symbolics.jacobian(C, sym_q)
Ct = expand_derivatives.(D.(C))

# 拘束力
Qc1 = Symbolics.jacobian(- Cq * sym_qdot, sym_q) * sym_qdot
Qc2 = expand_derivatives.(-2 * D.(Cq) * sym_qdot)
Qc3 = expand_derivatives.(- D.(Ct))
Qc = Qc1 + Qc2 + Qc3

# シンボリック表現を関数化
func_constraint = eval(build_function(C, sym_q)[1])
func_jacobian = eval(build_function(Cq, sym_q)[1])
func_Qc = eval(build_function(Qc, t, sym_q, sym_qdot)[1])


function G_bar(q::Vector)
    
    return SMatrix{3, 4}(2 * [
        -q[2] q[1] q[4] -q[3]
        -q[3] -q[4] q[1] q[2]
        -q[4] q[3] -q[2] q[1]
    ])

end

# システムの質量行列
function func_global_mass(q::AbstractVector)
    
    m_trans = diagm([m1, m1, m1])

    m_rot = transpose(G_bar(q[4:7])) * diagm([I1, I1, I1]) * G_bar(q[4:7])

    globalM = SMatrix{7, 7}(BlockDiagonal([m_trans, m_rot]))

    return globalM
end


"""
一般化外力
"""
function func_external_force()::SVector

    Q = SVector{7}([0, 0, -m1*g, 0, 0, 0, 0])
    
    return Q
end

function EOM(time, state)
    
    # 状態量をパースする
    q    = SVector{dim_coordinate}(state[1:dim_coordinate])
    qdot = SVector{dim_coordinate}(state[(dim_coordinate+1):end])

    # 一般化質量行列
    M = func_global_mass(q)

    # 一般化外力ベクトル
    Q = func_external_force()

    # 拘束条件式
    C = func_constraint(q)
    
    # ヤコビアン
    Cq = func_jacobian(q)
    
    Cdot = Cq * qdot

    Gm = func_Qc(0, q, qdot)

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
    Ts = 5e-5

    datanum = Integer(timelength / Ts + 1)
    times = 0.0:Ts:timelength
    states = [SVector{dim_coordinate * 2}(zeros(dim_coordinate * 2)) for _ in 1:datanum]
    states[1] = vcat([s1 * cos(pi/4), s1 * sin(pi/4), 0.0, 0.0, 0.0, 0.0, 0.0], zeros(dim_coordinate))
    accel = [SVector{dim_coordinate}(zeros(dim_coordinate)) for _ in 1:datanum]

    print("Running simulation ...")
    simtime = @elapsed for idx = 1:datanum-1

        # DAEから加速度を計算
        # accel[idx] = calc_acceleration(times[idx], states[idx])

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

(states, figures) = main();