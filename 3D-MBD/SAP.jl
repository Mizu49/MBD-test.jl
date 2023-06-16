using Revise, GLMakie, LinearAlgebra, StaticArrays, BlockDiagonals, GeometryBasics, Symbolics

include("RotationUtilities.jl")
using .RotationUtilities

# パラメータ設定
l1 = 1.0
s1 = l1/2
m1 = 0.5
I1 = m1 * l1^2 / 12

l2 = 1.0
s2 = l2/2
m2 = 0.5
I2 = m2 * l2^2 / 12

g = 9.81

# ボディの数と一般化座標の次元を定義
const num_bodies = 2
const dim_coordinate = 7
const dim_q = num_bodies * dim_coordinate

# シンボリック変数の定義
@variables t, sym_q[1:dim_q], sym_qdot[1:dim_q]
sym_q = collect(sym_q)
sym_qdot = collect(sym_qdot)

# オイラーパラメータ -> BRF
A_1 = transpose(quaternion2dcm(sym_q[4:7]))
A_2 = transpose(quaternion2dcm(sym_q[11:14]))

# Body 1の原点から拘束点への相対位置ベクトル
u_bar_1 = [0, -s1, 0]
u_bar_2 = [0, -s2, 0]

# 拘束条件のシンボリック表現
constraints_1 = [
    sym_q[1:3] + A_1 * u_bar_1 - (zeros(3)) # Revoluteジョイント位置に対する拘束
    transpose([1, 0, 0]) * (A_1 * [0, 1, 0]) # Revolute joint 1
    transpose([1, 0, 0]) * (A_1 * [0, 0, 1]) # Revolute joint 2
    transpose(sym_q[1:3]) * sym_q[1:3] - s1^2 # 絶対位置に対する拘束
]

constraints_2 = [
    sym_q[8:10] + A_2 * u_bar_2 - (sym_q[1:3] + A_1 * (-u_bar_1)) # Revoluteジョイント位置に対する拘束
    transpose([1, 0, 0]) * (A_2 * [0, 1, 0]) # Revolute joint 1
    transpose([1, 0, 0]) * (A_2 * [0, 0, 1]) # Revolute joint 2
    transpose(sym_q[8:10] - (sym_q[1:3] + A_1 * (-u_bar_1))) * (sym_q[8:10] - (sym_q[1:3] + A_1 * (-u_bar_1))) - s2^2 # 絶対位置に対する拘束（隣のボディとの結節点からの回転範囲内）
]


# 拘束条件をまとめる
C = vcat(constraints_1, constraints_2)

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

# システムの質量行列
function func_global_mass(q::AbstractVector)

    # body 1
    m_trans = diagm([m1, m1, m1])
    m_rot = transpose(G_bar(q[4:7])) * diagm([I1, I1, I1]) * G_bar(q[4:7])
    M_body1 = BlockDiagonal([m_trans, m_rot])

    # body 2
    m_trans = diagm([m2, m2, m2])
    m_rot = transpose(G_bar(q[11:14])) * diagm([I2, I2, I2]) * G_bar(q[11:14])
    M_body2 = BlockDiagonal([m_trans, m_rot])

    M_global = SMatrix{dim_q, dim_q}(BlockDiagonal([M_body1, M_body2]))

    return M_global
end


"""
一般化外力
"""
function func_external_force(time::Real, state::AbstractVector)::SVector

    # body 1
    body1_RPY = quaternion2euler(state[4:7])
    body1_force = [0.0, 0.0, 0.0]
    body1_torque = transpose(G(state[4:7])) * diagm([-0.25, -0.25, -0.25]) * body1_RPY
    # body1_torque = transpose(G_bar(state[4:7])) * [0.5, 0, 0]
    Q_body1 = vcat(body1_force, body1_torque)
    
    # body 2
    body2_RPY = quaternion2euler(state[11:14])
    body2_force = [0.0, 0.0, 0.0]
    # body2_torque = transpose(G(state[11:14])) * diagm([-5, -5, -5]) * body2_RPY
    body2_torque = transpose(G(state[11:14])) * zeros(3)
    Q_body2 = vcat(body2_force, body2_torque)

    Q = SVector{dim_q}(vcat(Q_body1, Q_body2))
    
    return Q
end

function EOM(time, state)
    
    # 状態量をパースする
    q    = SVector{dim_q}(state[1:dim_q])
    qdot = SVector{dim_q}(state[(dim_q+1):end])

    # 一般化質量行列
    M = func_global_mass(q)

    # 一般化外力ベクトル
    Q = func_external_force(time, state)

    # 拘束条件式
    C = func_constraint(q)
    
    # ヤコビアン
    Cq = func_jacobian(q)
    
    Cdot = Cq * qdot

    Gm = func_Qc(0, q, qdot)

    A = [
        M  transpose(Cq)
        Cq zeros(size(Cq, 1), size(Cq, 1))
    ]
    
    # バウムガルテの安定化法
    alpha = 10
    beta = 10
    Gm = Gm - 2 * alpha * Cdot - beta^2 * C

    # DAEの右辺を計算する
    RHS = vcat(Q, Gm)

    # 線形方程式を解いて一般化加速度とラグランジュ乗数を得る
    accel_lambda  = A \ RHS

    # 微分方程式の時間発展計算する 
    # qdot qddot
    differential = vcat(qdot, accel_lambda[1:dim_q])

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

    # 一般化座標＠グローバルの初期値計算
    init_body1_RPY = [pi/4, 0, 0] # 物体変換で規定
    init_body1_transposi = transpose(euler2dcm(init_body1_RPY)) * [0, 0.5, 0] # Localでの位置をGlobalでの位置に変換するときは座標変換行列の転置をつかう！
    init_body1_quaternion = euler2quaternion(init_body1_RPY) # 姿勢表現
    init_q1 = vcat(init_body1_transposi, init_body1_quaternion)

    init_body2_RPY = [pi/6, 0, 0]
    init_body2_transposi = init_body1_transposi + transpose(euler2dcm(init_body1_RPY)) * [0, 0.5, 0] + transpose(euler2dcm(init_body2_RPY)) * [0, 0.5, 0]
    init_body2_quaternion = euler2quaternion(init_body2_RPY)
    init_q2 = vcat(init_body2_transposi, init_body2_quaternion)

    init_q = vcat(init_q1, init_q2)

    # Static simulation with Newton-Raphson method
    print("Begin static analysis...")
    iter = 0
    isconv = false
    while (iter < 100) && isconv == false
        
        currentC  = func_constraint(init_q)
        currentCq = func_jacobian(init_q)

        dq = - pinv(currentCq) * currentC
        init_q = init_q + dq

        isconv_array = [false for _ in 1:size(currentC, 1)]
        for j = 1:(size(currentC, 1))
            if (abs(currentC[j, 1]) < 1e-10)
                isconv_array[j] = true 
            end
        end

        if all(isconv_array)
            is_conv = true
            println("complete!")
            break
        else
            iter = iter + 1
        end
    end

    println("initial coordinate: $init_q")

    # Dynamic simulation
    timelength = 10.0
    Ts = 1e-3

    datanum = Integer(timelength / Ts + 1)
    times = 0.0:Ts:timelength
    states = [SVector{dim_q * 2}(zeros(dim_q * 2)) for _ in 1:datanum]
    states[1] = vcat(init_q, zeros(dim_q))
    accel = [SVector{dim_q}(zeros(dim_q)) for _ in 1:datanum]

    print("Running dynamics simulation ...")
    simtime = @elapsed for idx = 1:datanum-1

        # DAEから加速度を計算
        # accel[idx] = calc_acceleration(times[idx], states[idx])

        # time evolution
        states[idx+1] = runge_kutta(times[idx], states[idx], Ts)

    end
    println("done!")

    # plot of x, y, z position for body 1
    fig1 = Figure()
    ax1 = Axis(
        fig1[1, 1],
        xlabel = "Time (s)",
        ylabel = "Position (m)"
    )

    x_position = getindex.(states, 1)
    y_position = getindex.(states, 2)
    z_position = getindex.(states, 3)

    lines!(ax1, times, x_position, label = "x-position")
    lines!(ax1, times, y_position, label = "y-position")
    lines!(ax1, times, z_position, label = "z-position")
    axislegend(ax1)

    # plot of x, y, z position for body 2
    fig2 = Figure()
    ax2 = Axis(
        fig2[1, 1],
        xlabel = "Time (s)",
        ylabel = "Position (m)"
    )

    x_position = getindex.(states, 8)
    y_position = getindex.(states, 9)
    z_position = getindex.(states, 10)

    lines!(ax2, times, x_position, label = "x-position")
    lines!(ax2, times, y_position, label = "y-position")
    lines!(ax2, times, z_position, label = "z-position")
    axislegend(ax2)

    figures = [fig1, fig2]

    return (times, states, figures)
end

(times, states, figures) = main();

Makie.inline!(true)
display(figures[1])
display(figures[2])
