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

# 拡大系の運動方程式に一般化座標を代入
substitute.(sysA, (Dict(q => rand(6)),))
