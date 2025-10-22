@kwdef struct FlipNSquare <: BH_Parameters{true,3,4,4}
    Lx::IndexType = 0
    Ly::IndexType = 0
    h0::f64 = 1.0
    h1::f64 = 1.0
    # bosonwl::Array{Wline,3} = [empty_wline(0, StateType(0)) for _ ∈ LinearIndices((Lx, Ly, 2))]
    buff::Wline = Element[]
    kink::Vector{Float64} = Float64[]
end

function wl_size(H::FlipNSquare)::NTuple{3,Int}
    return (Int(H.Lx), Int(H.Ly), 1)
end

function get_nbs(H::FlipNSquare, i::Integer)::NTuple{4,Int}
    @inbounds begin
        Lx, Ly, _ = wl_size(H)
        lid = LinearIndices((Lx, Ly))
        x0, y0 = Tuple(CartesianIndices(lid)[i])
        return (
            lid[mod1(x0 - 1, Lx), y0],  # left
            lid[mod1(x0 + 1, Lx), y0],  # right
            lid[x0, mod1(y0 - 1, Ly)],  # bottom
            lid[x0, mod1(y0 + 1, Ly)],  # top
        )
    end
end

function get_hps(H::FlipNSquare, i::Integer)::NTuple{4,Int}
    return get_nbs(H, i)
end

function diagE(H::FlipNSquare, i, ni::StateType, njs::NTuple{4,StateType})::f64
    return 0.0
end

function bond_weight(H::FlipNSquare, i::Integer, j::Integer)::f64
    return 0.0
end

function simple_measure_names(::Type{FlipNSquare})
    return (:E, :Sz, :Sz2, :K, :U, :μ, :V, :Kx, :Ky, :Wx², :Wy², :ρb, :Kb)
end

function simple_measure!(m,
    x::Wsheet{3},
    H::FlipNSquare
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
        NKb += length(l) - 1
    end
    ρb /= length(H.bosonwl)
    Kb = -NKb / x.β

    m.props = m.props .+ (E, abs2(E), N, abs2(N), K, U, μ, V, Kx, Ky, abs2(Wx), abs2(Wy), ρb, Kb)
    m.n_measure += 1
    return nothing
end

using Random
function update_boson!(H::Ham, x::Wsheet{N}, kb::Integer) where {Ham,N}
    i, j = bond_sites(H, kb)
    li = x[i]
    lb = H.bosonwl[kb]
    β = x.β
    buff = H.buff
    kink = H.kink
    empty!(kink)
    empty!(buff)
    for e ∈ li
        if e.j == j
            push!(kink, e.t)
        end
    end
    @assert issorted(kink)
    c = inv(β * H.Jb)
    μ = H.μb
    t = t_eps + c * randexp()
    n0 = lb[end].n_R
    for e ∈ lb
        while t < e.t - t_eps
            push!(buff, Element(t, IndexType(0), IndexType(0), n0, n0, I_))
            t += clamp(c * randexp(), t_eps, +Inf)
            if abs(t - e.t) < 2t_eps
                t += 5t_eps
            end
        end
        push!(buff, e)
        n0 = e.n_R
    end
    @assert issorted(buff)
    L = lastindex(buff) - 1
    modified = false
    if L ≤ 1 # try to flip the whole worldline
        @assert buff[1].op == I_
        if L == 1
            deleteat!(buff, 1)
        end
        @assert length(buff) == 1
        e = buff[1]
        n0 = e.n_R
        if ~isempty(kink)
            @assert n0 == StateType(1)
        end
        np = StateType(1 - n0)
        Pacc = exp(μ * β * (np - n0))
        if metro(Pacc) && isempty(kink)
            @assert e.t == 1.0
            buff[1] = Element(1.0, IndexType(0), IndexType(0), np, np, I_)
            modified = true
        end
    else
        for cur ∈ 1:L
            nxt = mod1(cur + 1, L)
            e_cur = buff[cur]
            e_nxt = buff[nxt]
            t1 = e_cur.t
            t2 = e_nxt.t
            Δt = mod(t2 - t1, 1.0)
            blk = if isempty(kink)
                false
            elseif t1 < t2
                searchsortedfirst(kink, t1) ≠ searchsortedfirst(kink, t2)
            else # if t2 < t1
                @assert cur == L
                first(kink) < t2 || last(kink) > t1
            end
            if blk
                continue
            end
            if ~(e_cur.n_R == e_nxt.n_L)
                println(cur => e_cur)
                println(nxt => e_nxt)
                println((cur, L, mod1(cur + 1, L)))
                println(buff)
                error("Wrong bosonwl")
            end
            n0 = e_cur.n_R
            np = StateType(1 - n0)
            dn = np - n0
            Pacc = exp(μ * β * Δt * dn)
            if metro(Pacc)
                buff[cur] = Element(e_cur.t, IndexType(0), IndexType(0),
                    e_cur.n_L, np, OpType(i8(e_cur.op) + dn)
                )
                buff[nxt] = Element(e_nxt.t, IndexType(0), IndexType(0),
                    np, e_nxt.n_R, OpType(i8(e_nxt.op) - dn)
                )
                if cur == L
                    buff[end] = Element(1.0, IndexType(0), IndexType(0),
                        np, np, I_
                    )
                end
                modified = true
            end
        end
    end
    if modified
        empty!(lb)
        for e ∈ buff
            if e.op ≠ I_
                push!(lb, e)
            end
        end
        push!(lb, last(buff))
        sizehint_wl!(lb)
    end
    return nothing
end

function check_config_bbq(x, H)
    for l ∈ x.wl
        for e ∈ l
            if 0 < e.i < e.j
                s = bond_state(H, e.i, e.j, e.t)
                @assert s == StateType(1)
            end
        end
    end
    return true
end