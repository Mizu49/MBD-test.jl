using LinearAlgebra, Plots, BenchmarkTools, StaticArrays

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

timelength = 500
Ts = 1e-3
datanum = Integer(timelength/Ts + 1)

initstate = [1, 0]

function simA(initstate, Ts, datanum)
    
    statesA = [zeros(2) for _ in 1:datanum]
    statesA[1] = initstate
    
    for idx = 1:datanum-1
        statesA[idx+1] = rk4(idx*Ts, statesA[idx], Ts)
    end

    return statesA
end

function simB(initstate, Ts, datanum)

    statesB = zeros(2, datanum)
    statesB[:, 1] = initstate
    for idx = 1:datanum-1 
        statesB[:, idx+1] = rk4(idx*Ts, statesB[:, idx], Ts)
    end
    
    return statesB
end

function simC(initstate, Ts, datanum)

    statesC = [SVector{2}(zeros(2)) for _ in 1:datanum]
    statesC[1] = initstate
    for idx = 1:datanum-1 
        statesC[idx+1] = rk4(idx*Ts, statesC[idx], Ts)
    end
    
    return statesC
end

function main()
    
    @time statesA = simA(initstate, Ts, datanum);
    @time statesB = simB(initstate, Ts, datanum);
    @time statesC = simC(initstate, Ts, datanum);

    return
end
