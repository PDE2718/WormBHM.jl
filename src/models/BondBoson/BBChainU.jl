@kwdef struct BBChainU <: BH_Parameters{true,2,2,2}
    nmax::StateType = 1
    L::IndexType = 0
    J::f64 = 1.0
    U::f64 = 0.0
    V::f64 = 0.0
    μ::f64 = 0.0
    μb::f64 = 1.0
    Ub::f64 = 1.0
    nBmax::Int = 1
    bosons::Array{Int,1} = zeros(Int, L)
end

function wl_size(H::BBChainU)::NTuple{2,Int}
    return (Int(H.L), 1)
end

function get_nbs(H::BBChainU, i::Integer)::NTuple{2,Int}
    return mod1.((i-1, i+1), Int(H.L))
end

function get_hps(H::BBChainU, i::Integer)::NTuple{2,Int}
    return get_nbs(H, i)
end

function diagE(H::BBChainU, i, ni::StateType, njs::NTuple{2,StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end

function bond_weight(H::BBChainU, i::Integer, j::Integer)::f64
    L = Int(H.L)
    j1, j2 = get_hps(H, i)
    return H.J * (
        if j == j1
            H.bosons[i]
        elseif j == j2
            H.bosons[j]
        else
            0
        end
    )
end

function bond_sites(H::BBChainU, kb::Integer)::NTuple{2,Int}
    return mod1.((kb,kb-1), Int(H.L))
end

function simple_measure_names(::Type{BBChainU})
    return (:E, :E², :N, :N², :K, :U, :μ, :V, :W², :Nb, :Nb²)
end

function simple_measure!(m,
    x::Wsheet{2},
    H::BBChainU
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
    Nb = sum(H.bosons)
    Nb2 = sum(abs2, H.bosons)
    # (:E, :E², :N, :N², :K, :U, :μ, :V, :W², :Nb, :Nb²)
    m.props = m.props .+ (E, abs2(E), N, abs2(N), K, U, μ, V, abs2(W), Nb, Nb2)
    m.n_measure += 1
    return nothing
end
