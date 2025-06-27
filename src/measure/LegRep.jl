using StaticArrays, Random
@kwdef mutable struct LegendreBin{N}
    const lb::Float64 = -1.0
    const ub::Float64 = +1.0
    counter::Int = 0
    Gl::SVector{N,Float64} = @SVector zeros(N)
end

@generated function collectPl_schmit(x::T, ::Val{N}, P_0::T=one(T)) where {N,T<:AbstractFloat}
    # @assert (N isa Integer) && (N > 0)
    return quote
        # P_0 = w # \sqrt 2*0+1
        P_1 = (√3) * x * P_0 # \sqrt 2*1+1
        Base.@nexprs $(N - 2) k -> begin
            P_{k + 1} = √(
                ((2k + 1) * (2k + 3)) // ((k + 1) * (k + 1))
            ) * x * P_k - √(
                (k * k * (2k + 3)) // ((k + 1) * (k + 1) * (2k-1))
            ) * P_{k - 1}
        end
        return SVector(Base.@ntuple $(N) k -> P_{k - 1})
    end
end

function accum_legbin!(B::LegendreBin{N}, τ::Float64, w::Float64 = 1.0) where {N}
    x = (2τ - B.lb - B.ub) / (B.ub - B.lb)
    @assert -1.0 ≤ x ≤ 1.0
    B.Gl += collectPl_schmit(x, Val{N}(), w)
    B.counter += 1
    return nothing
end

function interp(B::LegendreBin{N}, τ::Float64, lmax::Int = -1; norm::Bool = true, parity::Int=0) where {N}
    x = (2τ - B.lb - B.ub) / (B.ub - B.lb)
    @assert -1.0 ≤ x ≤ 1.0
    @assert lmax < N
    Gl = B.Gl
    Plx = collectPl_schmit(x, Val{N}(), 1.0)
    f = 0.
    for l ∈ 0:lmax
        if (parity==0) || iseven(l+parity)
            f += Plx[l+1] * Gl[l+1]
        end
    end
    if norm == true
        f /= B.counter
    end
    return f
end

import Base.merge
function merge(x::LegendreBin{N}, y::LegendreBin{N}) where {N}
    @assert x.lb == y.lb && x.ub == y.ub
    return LegendreBin{N}(
        x.lb, x.ub,
        x.counter + y.counter,
        x.Gl + y.Gl
    )
end