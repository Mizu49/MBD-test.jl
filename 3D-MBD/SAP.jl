using Revise, GLMakie, LinearAlgebra, StaticArrays, BlockDiagonals, GeometryBasics, Symbolics

include("rotation_utilities.jl")

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

dim_coordinate = 7

@variables t, sym_q[1:7], sym_qdot[1:7]
sym_q = collect(sym_q)
sym_qdot = collect(sym_qdot)

# オイラーパラメータ -> BRF
A_1 = [
    (sym_q[4]^2 - sym_q[5]^2 - sym_q[6]^2 + sym_q[7]^2) 2*(sym_q[4]*sym_q[5] + sym_q[6]*sym_q[7]) 2*(sym_q[4]*sym_q[6] - sym_q[5]*sym_q[7]);
    2*(sym_q[4]*sym_q[5] - sym_q[6]*sym_q[7]) (sym_q[5]^2 - sym_q[4]^2 - sym_q[6]^2 + sym_q[7]^2) 2*(sym_q[5]*sym_q[6] + sym_q[4]*sym_q[7]);
    2*(sym_q[4]*sym_q[6] + sym_q[5]*sym_q[7]) 2*(sym_q[5]*sym_q[6] - sym_q[4]*sym_q[7]) sym_q[6]^2 - sym_q[4]^2 - sym_q[5]^2 + sym_q[7]^2
]

# Body 1の原点から拘束点への相対位置ベクトル
u_bar_1 = [0, -s1, 0]

# 拘束条件のシンボリック表現
C = [
    sym_q[1:3] + A_1 * u_bar_1 - (zeros(3)) # 位置に対する拘束
    transpose([1, 0, 0]) * (A_1 * [0, 1, 0]) # Revolute joint 1
    transpose([1, 0, 0]) * (A_1 * [0, 0, 1]) # Revolute joint 2
    sym_q[2]^2 + sym_q[3]^2 - (s1)^2
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
        q[4] q[3] -q[2] -q[1]
        -q[3] q[4] q[1] -q[2]
        q[2] -q[1] q[4] -q[3]
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
function func_external_force(time::Real, state::AbstractVector)::SVector

    force = [0.1, 0.1, 0.1]
    torque = transpose(G_bar(state[4:7])) * [1, 0.05, 0.05]
    
    Q = SVector{7}(vcat(force, torque))
    
    return Q
end

function EOM(time, state)
    
    # 状態量をパースする
    q    = SVector{dim_coordinate}(state[1:dim_coordinate])
    qdot = SVector{dim_coordinate}(state[(dim_coordinate+1):end])

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
    differential = vcat(qdot, accel_lambda[1:dim_coordinate])

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

    # Static simulation
    println("Begin static analysis")

    dcm = C1(0.0)

    init_transposi = dcm * [0, 0.5, 0]
    init_eulerparam = dcm2quaternion(dcm)

    init_q = vcat(init_transposi, init_eulerparam)

    iter = 0
    isconv = false
    while (iter < 100) && isconv == false
        
        println(iter)

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
            println("finish!")
            break
        else
            iter = iter + 1
        end
    end

    println("init coordinate: $init_q")

    # Dynamic simulation


    timelength = 10.0
    Ts = 1e-3

    datanum = Integer(timelength / Ts + 1)
    times = 0.0:Ts:timelength
    states = [SVector{dim_coordinate * 2}(zeros(dim_coordinate * 2)) for _ in 1:datanum]
    states[1] = vcat(init_q, zeros(dim_coordinate))
    accel = [SVector{dim_coordinate}(zeros(dim_coordinate)) for _ in 1:datanum]

    print("Running dynamics simulation ...")
    simtime = @elapsed for idx = 1:datanum-1

        # DAEから加速度を計算
        # accel[idx] = calc_acceleration(times[idx], states[idx])

        # time evolution
        states[idx+1] = runge_kutta(times[idx], states[idx], Ts)

    end
    println("done!")

    # plot of x, y, z position
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

    figures = [fig1]

    return (times, states, figures)
end

(times, states, figures) = main();

Makie.inline!(true)
display(figures[1])
