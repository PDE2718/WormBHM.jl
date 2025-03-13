@kwdef struct BH_Square <: BH_Parameters{true, 3, 4, 4}
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 0 # 
    J::f64 = 1.0 # other parameters
    U::f64 = 0.0
    V::f64 = 0.0
    μ::f64 = 0.0
end
function wl_size(H::BH_Square)::NTuple{N_wldim(BH_Square),Int}
    return (Int(H.Lx),Int(H.Ly),1)
end
function get_nbs(H::BH_Square, i::Integer)::NTuple{N_nbs(BH_Square),Int}
    @inbounds begin
        Lx, Ly, _ = wl_size(H)
        lid = LinearIndices((Lx, Ly))
        x0, y0, = Tuple(CartesianIndices(lid)[i])
        return (
            lid[mod1(x0 - 1, Lx), y0],
            lid[mod1(x0 + 1, Lx), y0],
            lid[x0, mod1(y0 - 1, Ly)],
            lid[x0, mod1(y0 + 1, Ly)],
        )
    end
end
function get_hps(H::BH_Square, i::Integer)::NTuple{N_hps(BH_Square),Int}
    return get_nbs(H, i)
end
function diagE(H::BH_Square, i, ni::StateType, njs::NTuple{N_nbs(BH_Square), StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end
function bond_weight(H::BH_Square, i::Integer, j::Integer)::f64
    return (j ∈ get_hps(H, i)) ? H.J : 0.0
    # return H.J
end
function simple_measure_names(::Type{BH_Square})
    return (:E, :E², :N, :N², :K, :U, :μ, :V, :Kx, :Ky, :Wx², :Wy²)
end
function simple_measure!(m,
    x::Wsheet{N_wldim(BH_Square)},
    H::BH_Square
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
    # Wx, Wy, NKx, NKy = winding_number(x, H)
    Wx::Int = Wy::Int = 0
    NKx::Int = NKy::Int = 0
    for i ∈ eachindex(x)
        for e::Element ∈ x[i]::Wline
            if e.op == b_
                hps = get_hps(H, i)
                j = Int(e.j)
                if j == hps[1]
                    Wx -= 1
                    NKx += 1
                elseif j == hps[2]
                    Wx += 1
                    NKx += 1
                elseif j == hps[3]
                    Wy -= 1
                    NKy += 1
                elseif j == hps[4]
                    Wy += 1
                    NKy += 1
                else
                    error("illegal operator")
                end
            end
        end
    end
    @assert Wx % H.Lx == Wy % H.Ly == 0
    Wx = Wx ÷ H.Lx
    Wy = Wy ÷ H.Ly
    Kx = -NKx / x.β
    Ky = -NKy / x.β
    K = Kx + Ky
    E = U + μ + V + K
    # (:E, :E², :N, :N², :K, :U, :μ, :V, :Kx, :Ky, :Wx², :Wy²)
    m.props = m.props .+ (E, abs2(E), N, abs2(N), K, U, μ, V, Kx, Ky, abs2(Wx), abs2(Wy))
    m.n_measure += 1
    return nothing
end