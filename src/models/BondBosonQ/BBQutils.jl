using Random
function update_boson!(H::Ham, x::Wsheet{N}, kb::Integer) where {Ham, N}
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