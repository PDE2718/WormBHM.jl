@kwdef struct BBSquareU <: BH_Parameters{true,3,4,4}
    nmax::StateType = 1
    Lx::IndexType = 0
    Ly::IndexType = 0
    J::f64 = 1.0
    U::f64 = 0.0
    V::f64 = 0.0
    μ::f64 = 0.0
    μb::f64 = 1.0
    Ub::f64 = 1.0
    nBmax::Int = 1
    bosons::Array{Int,3} = zeros(Int, (Lx, Ly, 2))
end

function wl_size(H::BBSquareU)::NTuple{3,Int}
    return (Int(H.Lx), Int(H.Ly), 1)
end

function get_nbs(H::BBSquareU, i::Integer)::NTuple{4,Int}
    @inbounds begin
        Lx, Ly, _ = wl_size(H)
        lid = LinearIndices((Lx, Ly))
        x0, y0 = Tuple(CartesianIndices(lid)[i])
        return (
            lid[mod1(x0 - 1, Lx), y0],  # left
            lid[mod1(x0 + 1, Lx), y0],  # right
            lid[x0, mod1(y0 - 1, Ly)],  # bottom
            lid[x0, mod1(y0 + 1, Ly)],  # top
            # lid[x0, y0, mod1(z0 - 1, Lz)],  # back
            # lid[x0, y0, mod1(z0 + 1, Lz)],  # front
        )
    end
end

function get_hps(H::BBSquareU, i::Integer)::NTuple{4,Int}
    return get_nbs(H, i)
end

function diagE(H::BBSquareU, i, ni::StateType, njs::NTuple{4,StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end

function bond_weight(H::BBSquareU, i::Integer, j::Integer)::f64
    Lx, Ly = wl_size(H)
    hps = get_hps(H, i)
    x0, y0 = CartesianIndices((Lx, Ly))[i] |> Tuple
    return H.J * (
        if j == hps[1]  # left
            H.bosons[x0, y0, 1]
        elseif j == hps[2]  # right
            H.bosons[mod1(x0 + 1, Lx), y0, 1]
        elseif j == hps[3]  # bottom
            H.bosons[x0, y0, 2]
        elseif j == hps[4]  # top
            H.bosons[x0, mod1(y0 + 1, Ly), 2]
        else
            0
        end
    )
end

function bond_sites(H::BBSquareU, kb::Integer)::NTuple{2,Int}
    Lx, Ly = wl_size(H)
    x0, y0, sb = CartesianIndices((Lx, Ly, 2))[kb] |> Tuple
    i = LinearIndices((Lx,Ly))[x0, y0]
    return (i, get_hps(H, i)[2sb-1])
end

function simple_measure_names(::Type{BBSquareU})
    return (:E, :E², :N, :N², :K, :U, :μ, :V, :Kx, :Ky, :Wx², :Wy², :Nb, :Nb²)
end

function simple_measure!(m,
    x::Wsheet{3},
    H::BBSquareU
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
    Wx::Int = Wy::Int = 0
    NKx::Int = NKy::Int = 0
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
    Nb = sum(H.bosons)
    Nb2 = sum(abs2, H.bosons)
    m.props = m.props .+ (E, abs2(E), N, abs2(N), K, U, μ, V, Kx, Ky, abs2(Wx), abs2(Wy), Nb, Nb2)
    m.n_measure += 1
    return nothing
end

function BBDistBuffer(H::BBSquareU)
    return BBDistU[]
end

# function update_bosons!(H::BBSquareU, x::Wsheet{3}, fB::Vector{BBDistU})
#     B = H.bosons
#     a = x.β * 0.5 * H.Ub
#     b = -x.β * (0.5 * H.Ub + H.μb)
#     for kb ∈ eachindex(B)
#         i, j = bond_sites(H, kb)
#         Nk = count(e -> e.j == j, x[i])
#         np = rand_boson!(fB, Nk, a, b)
#         if 0 ≤ np ≤ H.nBmax
#             B[kb] = np
#         end
#     end
#     return nothing
# end

# function update_rand_boson!(H::BBSquareU, x::Wsheet{3}, fB::Vector{BBDistU})
#     B = H.bosons
#     a = x.β * 0.5 * H.Ub
#     b = -x.β * (0.5 * H.Ub + H.μb)
#     kij = rand(eachindex(B))
#     i, j = bond_sites(H, kij)
#     Nk = count(e -> e.j == j, x[i])
#     np = rand_boson!(fB, Nk, a, b)
#     if 0 ≤ np ≤ H.nBmax
#         B[kij] = np
#     end
#     return nothing
# end