mutable struct DensityMap{Nw}
    const ψ::Array{StateType,Nw}
    const ρ::Array{f64,Nw}
    n_measure::Int
end
function DensityMap(L::NTuple{Nw,<:Integer}) where {Nw}
    return DensityMap{Nw}(zeros(StateType, L), zeros(f64, L), 0)
end
DensityMap(x::Wsheet) = DensityMap(size(x))
function snapshot!(r::DensityMap{Nw}, x::Wsheet{Nw}, t_scaled::f64) where {Nw}
    # @assert size(r.ψ)  == size(r.ρ) == size(x.wl)
    @assert 0.0 ≤ t_scaled ≤ 1.0
    ψ = r.ψ
    wl = x.wl
    while t_scaled == 1.0
        t_scaled = rand()
    end
    if t_scaled == 0.0
        @inbounds for i ∈ eachindex(ψ, wl)
            ψ[i] = first(wl[i]).n_L
        end
    else
        @inbounds for i ∈ eachindex(ψ, wl)
            ψ[i] = slice_at(wl[i], t_scaled)
        end
    end
    r.ρ .+= ψ
    r.n_measure += 1
    return nothing
end
import Base.merge
function merge(ms::Array{DensityMap{Nw}}) where {Nw}
    ψ = similar(ms[1].ψ)
    ρ = sum(m.ρ for m ∈ ms)
    n_measure = sum(m.n_measure for m ∈ ms)
    return DensityMap{Nw}(
        ψ, ρ, n_measure
    )
end