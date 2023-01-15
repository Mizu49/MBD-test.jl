using Plots, LinearAlgebra, Symbolics, BlockDiagonals, StaticArrays

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
function M_system(element_masses::Vector{<:AbstractMatrix})
    # 要素数
    elemnum = size(element_masses, 1)
    
    # ブロック対角行列にする
    return SMatrix{3 * elemnum, 3 * elemnum}(BlockDiagonal(element_masses))
end

M1 = M_elem(m1, I1)
M2 = M_elem(m2, I2)
M = M_system([M1, M2])

@variables t, q[1:6]

# 拘束条件式
C = [
    q[1] - s1 * cos(q[3])
    q[2] - s1 * sin(q[3])
    q[1] + (l1 - s1) * cos(q[3]) - q[4] + s2 * cos(q[6])
    q[2] + (l1 - s1) * sin(q[3]) - q[5] + s2 * sin(q[6])
]
# 拘束条件式の次元
num_con = size(C, 1)

function calc_gamma(q, qdot)
    
    gamma = [
        -s1 * qdot[3]^2 * cos(q[3])
        -s1 * qdot[3]^2 * sin(q[3])
        (l1 - s1) * qdot[3]^2 * cos(q[3]) + s2 * qdot[6]^2 * cos(q[6])
        (l1 - s1) * qdot[3]^2 * sin(q[3]) + s2 * qdot[6]^2 * sin(q[6])
    ]

    return gamma
end


"""
    calc_Cq

シンボリック計算によりヤコビアンを計算する関数を返す
"""
function calc_Cq(C, q)

    # シンボリック計算によるヤコビアン計算
    Cq = Symbolics.jacobian(C, q)
    
    # 関数として返す
    expr = build_function(Cq, q)
    func = eval(expr[1])

    return func
end

# ヤコビアン
Cq = calc_Cq(C, q)

# 拡大系のシステムの運動方程式の構築
sysA = [
    M transpose(Cq(q))
    Cq(q) zeros(num_con, num_con)
]


function EOM(time, state)

    q = state[1:6]
    qdot = state[7:12]

    # 一般化外力ベクトル
    Q = [0, -m1*g, 0, 0, -m2*g, 0]

    # 拡大系の運動方程式に一般化座標を代入
    A = substitute.(sysA, (Dict(q => q),))

    Cdot = Cq(state[1:6]) * qdot
    
    gamma = calc_gamma(q, qdot)

    RHS = vcat(Q, gamma)

    accel_lambda  = A \ RHS

    # qdot qddot
    differential = vcat(qdot, accel_lambda[1:6])

    println(typeof(differential))

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
states = [zeros(6*2) for _ in 1:datanum]
states[1] = vcat([l1/2, 0.0, 0.0, l1 + l2/2, 0.0, 0.0], zeros(6))

for idx = 1:datanum-1
    
    # SVectorのVectorの要素ににVector{Num}を入れられないことがエラーの原因
    states[idx+1] = runge_kutta(times[idx], states[idx], Ts)

end

plt = plot()
plot!(plt, times, getindex.(states, 3), label = "phi1")
plot!(plt, times, getindex.(states, 6), label = "phi2")
xlabel!(plt, "Time (s)")
ylabel!(plt, "Angle (rad)")
display(plt)

