using LinearAlgebra, Plots

M = 10
C = 0
K = 100

sysA = [
    0    1
    -K/M -C/M
]

function EOM(time, state)
    return sysA * state
end

function rk4(time, state, Ts)

    k1 = EOM(time       , state            )
    k2 = EOM(time + Ts/2, state + Ts/2 * k1)
    k3 = EOM(time + Ts/2, state + Ts/2 * k2)
    k4 = EOM(time + Ts  , state + Ts   * k3)

    nextstate = state + Ts/6 * (k1 + 2 * k2 + 2 * k3 + k4)

    return nextstate
end

timelength = 10
Ts = 1e-2
datanum = Integer(timelength/Ts + 1)

initstate = [1, 0]
states = [zeros(2) for _ in 1:datanum]
states[1] = initstate

for idx = 1:datanum-1
    
    states[idx+1] = rk4(idx*Ts, states[idx], Ts)

end

plot(getindex.(states, 1))
plot!(getindex.(states, 2))
