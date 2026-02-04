mutable struct StructureFactorND{Ndim, Nsub, NSk}
    const ψs::NTuple{Nsub, Array{Float64, Ndim}}
    const ψk::NTuple{Nsub, Array{ComplexF64, Ndim}}
    const Sk::NTuple{NSk, Array{ComplexF64, Ndim}}
    const Sk_buf::Array{ComplexF64,Ndim}
    n_measure::Int
end
function rfft_buffer(L::NTuple{N,Int}) where {N}
    return Array{ComplexF64,N}(undef, (@set L[1] = 1 + L[1] >> 1))
end
function _sub_Sk_id(a::Integer, b::Integer, Nsub::Integer)
    return (2Nsub-a)*(a-1)÷2 + b
end
function _sub_Sk_eachindex(Nsub::Integer)
    return Base.OneTo(_sub_Sk_id(Nsub,Nsub,Nsub))
end
function empty!(S::StructureFactorND)
    for Ski ∈ S.Sk
        fill!(Ski, complex(0.0))
    end
    S.n_measure = 0
    return nothing
end
rfftplan(S::StructureFactorND) = plan_rfft(S.ψs |> first)

@generated function StructureFactorND(L::NTuple{Ndim, Integer}, ::Val{Nsub}) where {Ndim, Nsub}
    NSk = Nsub * (Nsub + 1) ÷ 2
    quote
        S = StructureFactorND{$(Ndim),$(Nsub),$(NSk)}(
            (@ntuple $(Nsub) i -> Array{Float64,$(Ndim)}(undef, L)),
            (@ntuple $(Nsub) i -> rfft_buffer(L)),
            (@ntuple $(NSk) i -> rfft_buffer(L)),
            # ntuple(i -> Array{Float64,$(Ndim)}(undef, L), Nsub),
            # ntuple(i -> rfft_buffer(L), Nsub),
            # ntuple(i -> rfft_buffer(L), NSk),
            rfft_buffer(L),
            0
        )
        empty!(S)
        return S
    end
end
function StructureFactorND(x::Wsheet{Nw}) where {Nw}
    L, Nsub = peel_last(size(x))
    return StructureFactorND(L, Val{Nsub}())
end
function StructureFactorND(x::Wsheet{Nw}, ::Val{Nsub}) where {Nw, Nsub}
    L, _Nsub = peel_last(size(x))
    @assert _Nsub == Nsub
    return StructureFactorND(L, Val{Nsub}())
end
function StructureFactorND(x::Array{T, Nw}) where {T, Nw}
    L, Nsub = peel_last(size(x))
    return StructureFactorND(L, Val{Nsub}())
end
function StructureFactorND(x::Array{T, Nw}, ::Val{Nsub}) where {T, Nw, Nsub}
    L, _Nsub = peel_last(size(x))
    @assert _Nsub == Nsub
    return StructureFactorND(L, Val{Nsub}())
end

function copy_states!(S::StructureFactorND, x::Array{T, Nw}) where {T, Nw}
    for (ψi, xi) ∈ zip(S.ψs, eachslice(x, dims = Nw))
        ψi .= xi
    end
    return nothing
end

function cal_Sk!(S::StructureFactorND{Ndim,Nsub,NSk}, P::FFTW.rFFTWPlan) where {Ndim, Nsub, NSk}
    for (ψk, ψ) ∈ zip(S.ψk, S.ψs)
        mul!(ψk, P, ψ)
    end
    i_Sk_count = 0
    for a ∈ 1:Nsub
        for b ∈ a:Nsub
            i_Sk_count += 1
            i_Sk = _sub_Sk_id(a, b, Nsub)
            @assert i_Sk_count == i_Sk
            map!(dot, S.Sk_buf, S.ψk[a], S.ψk[b])
            S.Sk[i_Sk_count] .+= S.Sk_buf
        end
    end
    @assert i_Sk_count == NSk
    S.n_measure += 1
    return nothing
end

@generated function measure_Sk!(S::StructureFactorND{Ndim,Nsub,NSk},
    ψsnap::Array{StateType, Nw}, P::FFTW.rFFTWPlan) where {Ndim,Nsub,NSk,Nw}
    @assert Nw == Ndim + 1
    return quote
        @assert peel_last(size(ψsnap)) == (size(first(S.ψs)), Nsub)
        @inbounds for i ∈ 1:Nsub
            S.ψs[i] .= view(ψsnap, (Base.@ntuple $(Ndim) i->Colon())..., i)
        end
        cal_Sk!(S, P)
        return nothing
    end
end
measure_Sk!(::Nothing, args...) = nothing

function Cab(S::StructureFactorND{Ndim,Nsub,NSk}, a::Integer, b::Integer) where {Ndim,Nsub,NSk}
    P = rfftplan(S)
    N = length(S.ψs |> first)
    @assert 0 < a ≤ Nsub && 0 < b ≤ Nsub
    if a ≤ b
        return inv(N * S.n_measure) * (P \ S.Sk[_sub_Sk_id(a,b,Nsub)])
    else
        return inv(N * S.n_measure) * (P \ conj.(S.Sk[_sub_Sk_id(b, a, Nsub)]))
    end
end
import Base: merge
function merge(ms::Array{StructureFactorND{Ndim, Nsub, NSk}}) where {Ndim, Nsub, NSk}
    m_merge = deepcopy(ms[1])
    for Sk ∈ m_merge.Sk
        Sk .= 0.0
    end
    for m ∈ ms
        for (Sk_bar, Sk) ∈ zip(m_merge.Sk, m.Sk)
            Sk_bar .+= Sk
        end
        m_merge.n_measure += m.n_measure
    end
    return m_merge
end