@kwdef struct BBChainQ <: BH_Parameters{true,2,2,2}
    nmax::StateType = 1
    L::IndexType = 0
    J::f64 = 1.0
    U::f64 = 0.0
    V::f64 = 0.0
    μ::f64 = 0.0
    μb::f64 = 0.0
    Jb::f64 = 1.0
    bosonwl::Array{Wline,2} = [empty_wline(0, StateType(0)) for _ ∈ LinearIndices((L, 1))]
    buff::Wline = Element[]
    kink::Vector{Float64} = Float64[]
end

function wl_size(H::BBChainQ)::NTuple{2,Int}
    return (Int(H.L), 1)
end

function get_nbs(H::BBChainQ, i::Integer)::NTuple{2,Int}
    return mod1.((i-1, i+1), Int(H.L))
end

function get_hps(H::BBChainQ, i::Integer)::NTuple{2,Int}
    return get_nbs(H, i)
end

function diagE(H::BBChainQ, i, ni::StateType, njs::NTuple{2,StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end

function bond_weight(H::BBChainQ, i::Integer, j::Integer)::f64
    return (j ∈ get_hps(H, i)) ? H.J : 0.0
end
function bond_state(H::BBChainQ, i::Integer, j::Integer, t::f64)::f64
    L = Int(H.L)
    j1, j2 = get_hps(H, i)
    lnb = if j == j1
            H.bosonwl[i]
        elseif j == j2
            H.bosonwl[j]
        else
            error("neighboring wl not found")
        end
    return slice_at(lnb, t)
end

function bond_sites(H::BBChainQ, kb::Integer)::NTuple{2,Int}
    return mod1.((kb,kb-1), Int(H.L))
end

function simple_measure_names(::Type{BBChainQ})
    return (:E, :E², :N, :N², :K, :U, :μ, :V, :W², :ρb, :Kb)
end

function simple_measure!(m,
    x::Wsheet{2},
    H::BBChainQ
)::Nothing
    U = μ = V = N = 0.0
    for (i, c) ∈ enumerate(CartesianIndices(x.wl))
        li::Wline = x[i]
        ni, ni² = integrated_density(Val(2), li)
        μ -= H.μ * ni
        U += 0.5 * H.U * (ni² - ni)
        for nb ∈ get_nbs(H, i)
            if nb < i
                V += H.V * integrated_densities(li, x[nb])
            end
        end
        N += li[end].n_L
    end
    @assert isapprox(μ, -H.μ * N)
    W = NK = 0
    for i ∈ eachindex(x)
        for e::Element ∈ x[i]::Wline
            if e.op == b_
                hps = get_hps(H, i)
                j = Int(e.j)
                if j == hps[1]  # left
                    W -= 1
                    NK += 1
                elseif j == hps[2]  # right
                    W += 1
                    NK += 1
                else
                    error("illegal operator")
                end
            end
        end
    end
    L = Int(H.L)
    @assert W % L == 0
    W = W ÷ L
    K = -NK / x.β
    E = U + μ + V + K

    ρb = 0.0
    NKb = 0
    for l ∈ H.bosonwl
        ρb += integrated_density(Val(1), l)[1]
        NKb += length(l) - 1
    end
    ρb /= length(H.bosonwl)
    Kb = -NKb / x.β
    m.props = m.props .+ (E, abs2(E), N, abs2(N), K, U, μ, V, abs2(W), ρb, Kb)
    m.n_measure += 1
    return nothing
end
