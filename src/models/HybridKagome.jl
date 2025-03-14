@kwdef struct HybridKagome <: BH_Parameters{true,3,4,4}
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 0 # So far only 2D is supported.
    J::f64 = 1.0
    ξ::f64 = 0.0
    V::f64 = 0.0
    μ::f64 = 0.0
end
function wl_size(H::HybridKagome)::NTuple{3,Int}
    return (Int(H.Lx), Int(H.Ly), 3)
end

function get_nbs(H::HybridKagome, i::Integer)::NTuple{4,Int}
    Lx = Int(H.Lx)
    Ly = Int(H.Ly)
    C = LinearIndices((Lx, Ly, 3))
    x0, y0, s0 = CartesianIndices((Lx, Ly, 3))[i] |> Tuple
    return if s0 == 1
        (
            C[x0, y0, 2],
            C[x0, y0, 3],
            C[mod1(x0 - 1, Lx), y0, 2],
            C[x0, mod1(y0 - 1, Ly), 3],
        )
    elseif s0 == 2
        (
            C[x0, y0, 3],
            C[x0, y0, 1],
            C[mod1(x0 + 1, Lx), mod1(y0 - 1, Ly), 3],
            C[mod1(x0 + 1, Lx), y0, 1],
        )
    else #if s0 == 3
        (
            C[x0, y0, 1],
            C[x0, y0, 2],
            C[x0, mod1(y0+1, Ly), 1],
            C[mod1(x0-1, Lx), mod1(y0 + 1, Ly), 2],
        )
    end
end
# The first two and the last two are for different plaquette.
get_hps(H::HybridKagome, i)::NTuple{4,Int} = get_nbs(H::HybridKagome, i)
function control_site(H::HybridKagome, i::Integer, j::Integer)::Int
    nbs = get_nbs(H, i)
    if j == nbs[1]
        return nbs[2]
    elseif j == nbs[2]
        return nbs[1]
    elseif j == nbs[3]
        return nbs[4]
    elseif j == nbs[4]
        return nbs[3]
    else
        error("wrong")
    end
end

function diagE(H::HybridKagome, i::Integer, ni::StateType, njs::NTuple{4,StateType})::f64
    return ni * (H.V * sum(njs) - H.μ)
end

function bond_weight(H::HybridKagome, i, j)::f64
    return (j ∈ get_hps(H, i)) ? H.J : 0.0
    # return H.J
end
function simple_measure_names(::Type{HybridKagome})
    return (:E, :E², :N, :N²,
        :N₁, :N₁², :N₂, :N₂², :N₃, :N₃²,
        :K, :μ, :V,
        :W₁², :W₂², :W₃², :W₁W₂, :W₂W₃, :W₃W₁,
    )
end
function simple_measure!(m, x::Wsheet{3}, H::HybridKagome)::Nothing
    μ = V = 0.0
    Ns = H.Lx * H.Ly
    N1 = N2 = N3 = 0
    W1 = W2 = W3 = 0
    K1 = K2 = K3 = 0
    for (i, c) ∈ enumerate(CartesianIndices(x.wl))
        li::Wline = x[i]
        ni, = integrated_density(Val(1), li)
        μ -= H.μ * ni
        for nb ∈ get_nbs(H, i)
            if nb < i
                V += H.V * integrated_densities(li, x[nb])
            end
        end
        if c[3] == 1
            N1 += li[1].n_L
            nbs = get_hps(H,i)
            for e ∈ li
                if e.op == b_
                    # @assert element_around(x.wl[control_site(H, e.i, e.j)], e.t, -1).n_R == 1
                    if e.j == nbs[1]
                        W1 += 2
                        W2 -= 1
                        W3 -= 1
                        K1 += 1
                    elseif e.j == nbs[2]
                        W3 -= 2
                        W1 += 1
                        W2 += 1
                        K3 += 1
                    elseif e.j == nbs[3]
                        W1 -= 2
                        W2 += 1
                        W3 += 1
                        K1 += 1
                    elseif e.j == nbs[4]
                        W3 += 2
                        W1 -= 1
                        W2 -= 1
                        K3 += 1
                    else
                        error("wrong kink")
                    end
                end
            end
        elseif c[3] == 2
            N2 += li[1].n_L
            nbs = get_hps(H,i)
            for e ∈ li
                if e.op == b_
                    # @assert element_around(x.wl[control_site(H, e.i, e.j)], e.t, -1).n_R == 1
                    if e.j == nbs[1]
                        W2 += 2
                        W3 -= 1
                        W1 -= 1
                        K2 += 1
                    elseif e.j == nbs[2]
                        W1 -= 2
                        W2 += 1
                        W3 += 1
                        K1 += 1
                    elseif e.j == nbs[3]
                        W2 -= 2
                        W3 += 1
                        W1 += 1
                        K2 += 1
                    elseif e.j == nbs[4]
                        W1 += 2
                        W2 -= 1
                        W3 -= 1
                        K1 += 1
                    else
                        error("wrong kink")
                    end
                end
            end
        elseif c[3] == 3
            N3 += li[1].n_L
            nbs = get_hps(H,i)
            for e ∈ li
                if e.op == b_
                    # @assert element_around(x.wl[control_site(H, e.i, e.j)], e.t, -1).n_R == 1
                    if e.j == nbs[1]
                        W3 += 2
                        W1 -= 1
                        W2 -= 1
                        K3 += 1
                    elseif e.j == nbs[2]
                        W2 -= 2
                        W3 += 1
                        W1 += 1
                        K2 += 1
                    elseif e.j == nbs[3]
                        W3 -= 2
                        W1 += 1
                        W2 += 1
                        K3 += 1
                    elseif e.j == nbs[4]
                        W2 += 2
                        W3 -= 1
                        W1 -= 1
                        K2 += 1
                    else
                        error("wrong kink")
                    end
                end
            end
        end
    end
    # @assert W1 % (H.Lx) == 0
    # W1 ÷= (H.Lx)
    # @assert W2 % (H.Ly) == 0
    # W2 ÷= (H.Ly)
    # @assert W3 % (H.Ly) == 0
    # W3 ÷= (H.Ly)
    # println("W1 = $(W1), W2 = $(W2), W3 = ($W3)")
    K = -(K1+K2+K3) / x.β
    N = N1 + N2 + N3
    E = μ + V + K
    # return (:E, :E², :N, :N²,
    #     :N₁, :N₁², :N₂, :N₂², :N₃, :N₃²,
    #     :K, :U, :μ, :V,
    #     :W₁², :W₂², :W₃², :W₁W₂, :W₂W₃, :W₃W₁,
    # )
    m.props = m.props .+ (
        E, abs2(E), N, abs2(N),
        N1, abs2(N1), N2, abs2(N2), N3, abs2(N3),
        K, μ, V,
        abs2(W1), abs2(W2), abs2(W3), W1*W2, W2*W3, W3*W1
    )
    m.n_measure += 1
    return nothing
end
