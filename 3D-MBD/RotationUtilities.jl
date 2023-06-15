module RotationUtilities

using LinearAlgebra, Symbolics, StaticArrays

export C1, C2, C3, G, G_bar, dcm2quaternion, dcm2euler, euler2dcm, euler2quaternion, quaternion2dcm, quaternion2euler

"""
    C1(theta::Real)::SMatrix

Rotational matrix for 1-axis
"""
@inline function C1(theta::Real)::SMatrix
    return SMatrix{3, 3, <:Real}([
        1 0 0
        0 cos(theta) -sin(theta)
        0 sin(theta) cos(theta)
    ])
end

@inline function C1(theta::Num)::Matrix{Num}
    return [
        1 0 0
        0 cos(theta) -sin(theta)
        0 sin(theta) cos(theta)
    ]
end


"""
    C2(theta::Real)::SMatrix

Rotational matrix for 2-axis
"""
@inline function C2(theta::Real)::SMatrix
    return SMatrix{3, 3, <:Real}([
        cos(theta) 0 sin(theta)
        0 1 0
        -sin(theta) 0 cos(theta)
    ])
end

@inline function C2(theta::Num)::Matrix{Num}
    return [
        cos(theta) 0 sin(theta)
        0 1 0
        -sin(theta) 0 cos(theta)
    ]
end

"""
    C3(theta::Real)::SMatrix

Rotational matrix for 3-axis
"""
@inline function C3(theta::Real)::SMatrix
    return SMatrix{3, 3, <:Real}([
        cos(theta) -sin(theta) 0
        sin(theta) cos(theta) 0
        0 0 1
    ])
end

@inline function C3(theta::Num)::Matrix{Num}
    return [
        cos(theta) -sin(theta) 0
        sin(theta) cos(theta) 0
        0 0 1
    ]
end

@inline function Base.:~(x::AbstractVector)
    return SMatrix{3, 3, <:Real}([
        0 -x[3] x[2]
        x[3] 0 -x[1]
        -x[2] x[1] 0
    ])
end

@inline function E(q::Vector{<:Real})
    return SMatrix{3, 4}(hcat(q[4] * I + ~(q[1:3]), -q[1:3]))
end

@inline function E_bar(q::Vector{<:Real})
    return SMatrix{3, 4}(hcat(q[4] * I - ~(q[1:3]), -q[1:3]))
end

@inline function G(q::Vector)
    return 2 * E(q)
end

@inline function G_bar(q::Vector)
    return 2 * E_bar(q)
end


"""
    dcm2quaternion(dcm::Matrix{Real})::Vector{Real}

calculate quaternion from direction cosine matrix (DCM) `dcm`
"""
@inline function dcm2quaternion(dcm::SMatrix{3, 3, Float64})::SVector{4, Real}

    _checkdcm(dcm)

    q = [
        sqrt(1 + dcm[1,1] - dcm[2,2] - dcm[3,3])/2,
        sqrt(1 - dcm[1,1] + dcm[2,2] - dcm[3,3])/2,
        sqrt(1 - dcm[1,1] - dcm[2,2] + dcm[3,3])/2,
        sqrt(1 + dcm[1,1] + dcm[2,2] + dcm[3,3])/2
    ]

    (maxvalue, maxindex) = findmax(q)

    if maxindex == 1
        q[2] = 0.25/q[1] * (dcm[1,2] + dcm[2,1])
        q[3] = 0.25/q[1] * (dcm[1,3] + dcm[3,1])
        q[4] = 0.25/q[1] * (dcm[2,3] - dcm[3,2])
    elseif maxindex == 2
        q[1] = 0.25/q[2] * (dcm[1,2] + dcm[2,1])
        q[3] = 0.25/q[2] * (dcm[3,2] + dcm[2,3])
        q[4] = 0.25/q[2] * (dcm[3,1] - dcm[1,3])
    elseif maxindex == 3
        q[1] = 0.25/q[3] * (dcm[3,1] + dcm[1,3])
        q[2] = 0.25/q[3] * (dcm[3,2] + dcm[2,3])
        q[4] = 0.25/q[3] * (dcm[1,2] - dcm[2,1])
    elseif maxindex == 4
        q[1] = 0.25/q[4] * (dcm[2,3] - dcm[3,2])
        q[2] = 0.25/q[4] * (dcm[3,1] - dcm[1,3])
        q[3] = 0.25/q[4] * (dcm[1,2] - dcm[2,1])
    else
        error("`maxindex` is illegal")
    end

    return SVector{4}(q)
end

@inline function dcm2quaternion(dcm::Matrix{Num})::Vector{Num}

    error("This function is still not implemented!")

    return 
end

"""
    function eular2dcm(euler::Union{SVector{3, <:Real}, Vector{<:Real}})::SMatrix{3, 3, Float64}

calculate direction cosine matrix from the vector of z-y-x eular angles.

## Argument

* `euler::Union{SVector{3, <:Real}, Vector{<:Real}}`: each element represents the rotation with z, y, x axis, respectively
"""
@inline function euler2dcm(euler::Union{SVector{3, <:Real}, Vector{<:Real}})::SMatrix{3, 3, Float64}
    return C3(euler[3]) * C2(euler[2]) * C1(euler[1])
end

"""
    quaternion2dcm(q::Union{Vector{<:Real}, SVector{4, <:Real}})

calculates direction cosine matrix from quaternion
"""
@inline function quaternion2dcm(q::Union{Vector{<:Real}, SVector{4, <:Real}})::SMatrix{3, 3, Float64}
    q2 = q.^2;
    dcm = SMatrix{3, 3}([
        (q2[1] - q2[2] - q2[3] + q2[4]) 2*(q[1]*q[2] + q[3]*q[4]) 2*(q[1]*q[3] - q[2]*q[4]);
        2*(q[1]*q[2] - q[3]*q[4]) (q2[2] - q2[1] - q2[3] + q2[4]) 2*(q[2]*q[3] + q[1]*q[4]);
        2*(q[1]*q[3] + q[2]*q[4]) 2*(q[2]*q[3] - q[1]*q[4]) q2[3] - q2[1] - q2[2] + q2[4]
    ])
    return dcm
end

@inline function quaternion2dcm(q::Vector{Num})::Matrix{Num}
    q2 = q.^2;
    dcm = [
        (q2[1] - q2[2] - q2[3] + q2[4]) 2*(q[1]*q[2] + q[3]*q[4]) 2*(q[1]*q[3] - q[2]*q[4]);
        2*(q[1]*q[2] - q[3]*q[4]) (q2[2] - q2[1] - q2[3] + q2[4]) 2*(q[2]*q[3] + q[1]*q[4]);
        2*(q[1]*q[3] + q[2]*q[4]) 2*(q[2]*q[3] - q[1]*q[4]) q2[3] - q2[1] - q2[2] + q2[4]
    ]
    return dcm
end


"""
    dcm2euler(dcm::Union{SMatrix{3, 3, <:Real}, Matrix{<:Real}})::SVector{3, <:Real}

calculates z-y-x euler rotation angle from direction cosine matrix
"""
@inline function dcm2euler(dcm::Union{SMatrix{3, 3, <:Real}, Matrix{<:Real}})::SVector{3, <:Real}
    _checkdcm(dcm)
    # 3-2-1 euler angle (roll-pitch-yaw)
    euler = SVector{3}(-[
        atan(dcm[2,3], dcm[3,3]),
        atan(-dcm[1,3], sqrt(dcm[2,3]^2 + dcm[3,3]^2)),
        atan(dcm[1,2], dcm[1,1])
    ])
    return euler
end

@inline function dcm2euler(dcm::Matrix{Num})::Vector{Num}
    _checkdcm(dcm)
    # 3-2-1 euler angle (roll-pitch-yaw)
    euler = -[
        atan(dcm[2,3], dcm[3,3]),
        atan(-dcm[1,3], sqrt(dcm[2,3]^2 + dcm[3,3]^2)),
        atan(dcm[1,2], dcm[1,1])
    ]
    return euler
end

"""
    quaternion2euler(quaternion::Union{Vector{<:Real}, SVector{4, <:Real}})::SVector{3, <:Real}

calculates z-y-x euler rotation angle from quaternion
"""
@inline function quaternion2euler(quaternion::Union{Vector{<:Real}, SVector{4, <:Real}})::SVector{3, <:Real}

    # use DCM for the calculation
    euler = dcm2euler(quaternion2dcm(quaternion))

    return euler
end

@inline function quaternion2euler(quaternion::Vector{Num})::Vector{Num}

    # use DCM for the calculation
    euler = dcm2euler(quaternion2dcm(quaternion))

    return euler
end

"""
    euler2quaternion(euler::Union{SVector{3, <:Real}, Vector{<:Real}})::SVector{4, Real}

calculates quaternion from z-y-x euler rotation angle
"""
@inline function euler2quaternion(euler::Union{SVector{3, <:Real}, Vector{<:Real}})::SVector{4, Real}

    # use DCM for the calculation
    quaternion = dcm2quaternion(euler2dcm(euler))

    return quaternion
end

@inline function euler2quaternion(euler::Vector{Num})::Vector{Num}

    # use DCM for the calculation
    quaternion = dcm2quaternion(euler2dcm(euler))

    return quaternion
end

@inline function _checkdcm(dcm::Union{SMatrix{3, 3, <:Real}, Matrix{<:Real}})
    if size(dcm) != (3, 3)
        throw(ArgumentError("`dcm` should be `3x3` matrix"))
    end
    return
end

"""
    Base.Math.deg2rad(rotationangle::Union{SVector{3, <:Real}, Vector{<:Real})

Convert rotation angle vector in degrees to radians
"""
function Base.Math.deg2rad(rotationangle::Union{SVector{3, <:Real}, Vector{<:Real}})

    newangle = SVector{3}([
        deg2rad(rotationangle[1])
        deg2rad(rotationangle[2])
        deg2rad(rotationangle[3])
    ])

    return newangle
end

"""
    Base.Math.rad2deg(rotationangle::Union{SVector{3, <:Real}, Vector{<:Real})

Convert rotation angle vector in radians to degrees
"""
function Base.Math.rad2deg(rotationangle::Union{SVector{3, <:Real}, Vector{<:Real}})

    newangle = SVector{3}([
        rad2deg(rotationangle[1])
        rad2deg(rotationangle[2])
        rad2deg(rotationangle[3])
    ])

    return newangle
end

end
