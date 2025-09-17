@kwdef struct BBSquareQ <: BH_Parameters{true,3,4,4}
    nmax::StateType = 1
    Lx::IndexType = 0
    Ly::IndexType = 0
    J::f64 = 1.0
    U::f64 = 0.0
    V::f64 = 0.0
    μ::f64 = 0.0
    μb::f64 = 0.0
    Jb::f64 = 1.0
    bosonwl::Array{Wline, 3} = [empty_wline(0, StateType(0)) for _ ∈ LinearIndices((Lx, Ly, 2))]
    buff::Wline = Element[]
    kink::Vector{Float64} = Float64[]
end

function wl_size(H::BBSquareQ)::NTuple{3,Int}
    return (Int(H.Lx), Int(H.Ly), 1)
end

function get_nbs(H::BBSquareQ, i::Integer)::NTuple{4,Int}
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

function get_hps(H::BBSquareQ, i::Integer)::NTuple{4,Int}
    return get_nbs(H, i)
end

function diagE(H::BBSquareQ, i, ni::StateType, njs::NTuple{4,StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end

function bond_weight(H::BBSquareQ, i::Integer, j::Integer)::f64
    return (j ∈ get_hps(H, i)) ? H.J : 0.0
    # return H.J
end
function bond_state(H::BBSquareQ, i::Integer, j::Integer, t::f64)::f64
    Lx, Ly = wl_size(H)
    hps = get_hps(H, i)
    x0, y0 = CartesianIndices((Lx, Ly))[i] |> Tuple
    lnb = if j == hps[1]  # left
            H.bosonwl[x0, y0, 1]
        elseif j == hps[2]  # right
            H.bosonwl[mod1(x0 + 1, Lx), y0, 1]
        elseif j == hps[3]  # bottom
            H.bosonwl[x0, y0, 2]
        elseif j == hps[4]  # top
            H.bosonwl[x0, mod1(y0 + 1, Ly), 2]
        else
            error("neighboring wl not found")
        end
    return slice_at(lnb, t)
end

function bond_sites(H::BBSquareQ, kb::Integer)::NTuple{2,Int}
    Lx, Ly = wl_size(H)
    x0, y0, sb = CartesianIndices((Lx, Ly, 2))[kb] |> Tuple
    i = LinearIndices((Lx,Ly))[x0, y0]
    return (i, get_hps(H, i)[2sb-1])
end

function simple_measure_names(::Type{BBSquareQ})
    return (:E, :E², :N, :N², :K, :U, :μ, :V, :Kx, :Ky, :Wx², :Wy², :ρb, :Kb)
end

function simple_measure!(m,
    x::Wsheet{3},
    H::BBSquareQ
)::Nothing
    U = μ = V = N = 0.0
    for (i, c) ∈ enumerate(CartesianIndices(x.wl))
        li::Wline = x[i]
        ni, ni² = integrated_density(Val(2), li)
        μ -= H.μ * ni
        U += 0.5 * H.U * (ni² - ni)
        if H.V ≠ 0.0
            for nb ∈ get_nbs(H, i)
                if nb < i
                    V += H.V * integrated_densities(li, x[nb])
                end
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

    ρb = 0.0
    NKb = 0
    for l ∈ H.bosonwl
        ρb += integrated_density(Val(1), l)[1]
        NKb += length(l)-1
    end
    ρb /= length(H.bosonwl)
    Kb = - NKb / x.β

    m.props = m.props .+ (E, abs2(E), N, abs2(N), K, U, μ, V, Kx, Ky, abs2(Wx), abs2(Wy), ρb, Kb)
    m.n_measure += 1
    return nothing
end