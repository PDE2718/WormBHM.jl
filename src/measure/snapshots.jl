mutable struct FixedSizeRecorder{T}
    const L::Int
    t::Int
    Δ::Int
    snaps::Vector{T}
end
function FixedSizeRecorder(::Type{T}, L::Integer) where {T}
    @assert L > 0 && iseven(L)
    return FixedSizeRecorder{T}(Int(L), 0, 1, T[])
end
function record!(x::FixedSizeRecorder{T}, v::T) where {T}
    x.t += 1
    if x.t % x.Δ == 0
        if length(x.snaps) < x.L
            push!(x.snaps, deepcopy(v))
        else
            x.Δ *= 2
            deleteat!(x.snaps, 1:2:length(x.snaps))
        end
    end
    nothing
end
record!(snapshots, density_map.ψ)

function SnapShots(x::Wsheet{Nw}, L::Integer) where {Nw}
    return FixedSizeRecorder(Array{StateType, Nw}, L)
end

mutable struct CoarseGrainedRecorder{T}
    l::Vector{Vector{T}}
end
function CoarseGrainedRecorder(T::DataType)
    return CoarseGrainedRecorder([T[]])
end
import Base: <<, empty!
function <<(r::CoarseGrainedRecorder{T}, val::T) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}
    l = r.l
    push!(l[1], val)
    @inbounds for i ∈ 1:(length(l)-1)
        li = l[i]
        lj = l[i+1]
        L = length(li)
        if L |> iseven
            push!(lj, (li[L-1] + li[L]) / 2)
        else
            break
        end
    end
    if length(l[end]) == 2
        push!(l, T[(l[end][1]+l[end][2])/2])
    end
    return r
end