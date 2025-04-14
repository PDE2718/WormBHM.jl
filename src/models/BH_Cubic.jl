@kwdef struct BH_Cubic <: BH_Parameters{true, 4, 6, 6}
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 0 # 
    Lz::IndexType = 0 #
    J::f64 = 1.0 # other parameters
    U::f64 = 0.0
    V::f64 = 0.0
    μ::f64 = 0.0
end

function wl_size(H::BH_Cubic)::NTuple{4,Int}
    return (Int(H.Lx), Int(H.Ly), Int(H.Lz), 1)
end

function get_nbs(H::BH_Cubic, i::Integer)::NTuple{6,Int}
    @inbounds begin
        Lx, Ly, Lz, _ = wl_size(H)
        lid = LinearIndices((Lx, Ly, Lz))
        x0, y0, z0 = Tuple(CartesianIndices(lid)[i])
        return (
            lid[mod1(x0 - 1, Lx), y0, z0],  # left
            lid[mod1(x0 + 1, Lx), y0, z0],  # right
            lid[x0, mod1(y0 - 1, Ly), z0],  # bottom
            lid[x0, mod1(y0 + 1, Ly), z0],  # top
            lid[x0, y0, mod1(z0 - 1, Lz)],  # back
            lid[x0, y0, mod1(z0 + 1, Lz)],  # front
        )
    end
end

function get_hps(H::BH_Cubic, i::Integer)::NTuple{6,Int}
    return get_nbs(H, i)
end

function diagE(H::BH_Cubic, i, ni::StateType, njs::NTuple{6,StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end

function bond_weight(H::BH_Cubic, i::Integer, j::Integer)::f64
    return (j ∈ get_hps(H, i)) ? H.J : 0.0
    # return H.J
end

function simple_measure_names(::Type{BH_Cubic})
    return (:E, :E², :N, :N², :K, :U, :μ, :V, :Kx, :Ky, :Kz, :Wx², :Wy², :Wz²)
end

function simple_measure!(m,
    x::Wsheet{4},
    H::BH_Cubic
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
    Wx::Int = Wy::Int = Wz::Int = 0
    NKx::Int = NKy::Int = NKz::Int = 0
    for i ∈ eachindex(x)
        for e::Element ∈ x[i]::Wline
            if e.op == b_
                hps = get_hps(H, i)
                j = Int(e.j)
                if j == hps[1]  # left
                    Wx -= 1
                    NKx += 1
                elseif j == hps[2]  # right
                    Wx += 1
                    NKx += 1
                elseif j == hps[3]  # bottom
                    Wy -= 1
                    NKy += 1
                elseif j == hps[4]  # top
                    Wy += 1
                    NKy += 1
                elseif j == hps[5]  # back
                    Wz -= 1
                    NKz += 1
                elseif j == hps[6]  # front
                    Wz += 1
                    NKz += 1
                else
                    error("illegal operator")
                end
            end
        end
    end
    @assert Wx % H.Lx == Wy % H.Ly == Wz % H.Lz == 0
    Wx = Wx ÷ H.Lx
    Wy = Wy ÷ H.Ly
    Wz = Wz ÷ H.Lz
    Kx = -NKx / x.β
    Ky = -NKy / x.β
    Kz = -NKz / x.β
    K = Kx + Ky + Kz
    E = U + μ + V + K
    # return (:E, :E², :N, :N², :K, :U, :μ, :V, :Kx, :Ky, :Kz, :Wx², :Wy², :Wz²)
    m.props = m.props .+ (E, abs2(E), N, abs2(N), K, U, μ, V, Kx, Ky, Kz, abs2(Wx), abs2(Wy), abs2(Wz))
    m.n_measure += 1
    return nothing
end