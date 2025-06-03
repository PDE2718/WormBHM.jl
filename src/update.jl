@generated function worm_cycle!(
        x::Wsheet{Nw},
        H::Ham,
        Cw::f64, Eoff::f64,
        AP_table::NTuple{4, f64},
        n_cycle::Int,
        G::T_G = nothing,
        fB::T_FB = nothing,
    )::Int where {
        Nw, Ham<:BH_Parameters,
        T_G<:Union{GreenFuncBin, Nothing},
        T_FB<:Union{BondSampler, Nothing}
    }
    
    @assert Nw == N_wldim(Ham)
    znbs::Int = N_nbs(Ham)
    zhps::Int = N_hps(Ham)
    mes_green::Bool = ~(T_G == Nothing)
    is_KCM_Kagome::Bool = (Ham == HybridKagome)
    is_Bond_Boson::Bool = (:bosons ∈ fieldnames(Ham))
    # println("compiling update for $(Ham) (Ndim = $(Ndim), znbs = $(znbs)), zhps = $(zhps), GF = $(mes_green)")
    quote
        $(
            if mes_green
                quote
                    G.insertion_trial += n_cycle
                end 
            end
        )
        cycle_size::Int = 0
        lattice = Base.OneTo(IndexType(lastindex(x.wl)))
        F_move, F_ins, F_del, F_glue = AP_table
        $(
            if is_Bond_Boson
                quote
                    @assert F_glue < 1.0
                end
            else
                quote
                    @assert F_glue == 1.0
                end
            end
        )
        P_del2ins = (F_del - F_ins)/(F_ins - F_move)
        P_glue = F_glue - F_del

        @inbounds for _ ∈ 1:n_cycle
            #############################################################################
            ############################# [INSERT_WORM] #################################
            #############################################################################
            # relative direction from b to b̂ / other direction choice
            δ::StateType = D::StateType = randsign()
            loc::WormLocation = _at_stop
            i::IndexType = rand(lattice)
            li::Wline = x[i]
            t::f64 = rand()
            head_id::Int = vindex(li, t)
            n0::StateType = li[head_id].n_L
            np::StateType = n0 - D
            # P_acc::f64 = 2 * Cw * max(np, n0)
            P_acc::f64 = 2 * Cw * max(np, n0) * P_glue
            # case insertion failed!
            if np < StateType(0) || np > H.nmax || t < t_eps || t > (1 - t_eps) || !(metro(P_acc)) || close_to_any(li, t)
                continue
            end
            hps::NTuple{$zhps,Int} = get_hps(H, i)::NTuple{$zhps,Int}
            @nexprs $zhps k -> begin
                if hps[k] > 0 && close_to_any(x[hps[k]], t)
                    continue
                end
            end # otherwise accepted
            nbs::NTuple{$zhps,Int} = get_nbs(H, i)::NTuple{$znbs,Int}
            tail::Element = Element(t, i, IndexType(0),
                D == i8(+1) ? n0 : np,
                D == i8(+1) ? np : n0,
                b_
            )
            head::Element = Element(nextfloat(t, D), i, StateType(0), tail.n_R, tail.n_L, b̂_)
            $(
                if is_KCM_Kagome
                    quote
                        kcm_nb_kink::Element = Element()
                    end
                end
            )
            insert!(li, head_id, tail)
            insert!(li, head_id + (D == i8(+1)), head)
            if rand(Bool) # swap head and tail
                tail, head = (head, tail)
                δ = -δ
            end
            head_id = vindex(li, head.t)
            cycle_size += 1
            sizehint_wl!(li)
            
            @label CYCLE_START🔁

            $(
                if mes_green
                    quote
                        if loc ∈ (_at_free, _at_green, _at_stop)
                            accumulate!(G, x, tail, head, loc)
                        end
                    end
                end
            )
            dice = rand()
            if dice < F_move     # [MOVE_WORM]
                # println("   moving...")
                D = randsign()
                if D * δ == i8(-1)
                    if loc == _at_stop || loc == _at_kink
                        @goto CYCLE_START🔁
                    elseif loc == _at_nbkink
                        $(
                            if is_KCM_Kagome
                                quote
                                    # @assert nbs == get_nbs(H, head.i)
                                    # @assert kcm_nb_kink.i ∈ nbs
                                    # @assert kcm_nb_kink.t == nextfloat(head.t, D)
                                    if kcm_nb_kink.j ∈ nbs
                                        # @assert kcm_nb_kink.n_L ≠ kcm_nb_kink.n_R
                                        n_k = D == i8(+1) ? head.n_R : head.n_L
                                        P_pass = n_k == StateType(1) ? H.ξ : inv(H.ξ)
                                        if metro(P_pass)
                                            head <<= nextfloat(head.t, 2D)
                                            li[head_id] = head
                                            δ = -δ
                                        end
                                        @goto CYCLE_START🔁
                                    end
                                end
                            end
                        )
                        head <<= nextfloat(head.t, 2D)
                        li[head_id] = head
                        δ = -δ
                    elseif loc == _at_green
                        head <<= nextfloat(head.t, 2D)
                        li[head_id] = head
                        δ = -δ
                    elseif loc == _at_dummy # pass dummy
                        dummy = li[end]
                        if δ == i8(+1) && D == i8(-1)
                            # @assert head == li[1] && head.t == nextfloat(0.)
                            head <<= prevfloat(1.0)
                            dummy = dummy_element(dummy.i, head.n_R)
                            li[end] = head
                            push!(li, dummy)
                            popfirst!(li)
                            head_id = lastindex(li) - 1
                        else  # if δ == -1 && D == +1
                            head <<= nextfloat(0.0)
                            dummy = dummy_element(dummy.i, head.n_L)
                            pushfirst!(li, head)
                            pop!(li)
                            li[end] = dummy
                            head_id = 1
                        end
                        δ = -δ
                    end
                end
                t = head.t
                v_near::Element = li[mod1(head_id + D, length(li))]
                # nbs::NTuple{$znbs,Int} = get_nbs(H, i)::NTuple{$znbs,Int}
                @nexprs $znbs k -> begin
                    assoc_k = if nbs[k] > 0
                        element_around(x[nbs[k]], t, D)
                    else
                        Element()
                    end
                end
                nb_states = if D == StateType(+1)
                    @ntuple $znbs k -> assoc_k.n_L
                else
                    @ntuple $znbs k -> assoc_k.n_R
                end
                λ₊::f64 = λ₋::f64 = ΔE::f64 = Eoff
                Eb::f64 = diagE(H, i, head.n_L, nb_states)
                Ef::f64 = diagE(H, i, head.n_R, nb_states)
                # println("Eb,Ef = $((Eb,Ef))")
                if D == StateType(-1)
                    ΔE = Ef
                    Ef = Eb
                    Eb = ΔE
                end
                ΔE = abs(Ef - Eb)
                Δt = randexp()
                ##################################### try to move
                if D == StateType(+1) # move forward -->
                    @nexprs $znbs k -> begin
                        if assoc_k.i > IndexType(0) && assoc_k.t < v_near.t
                            v_near = assoc_k
                        end
                    end
                    # however, the head cannot pass the tail
                    if head.t < tail.t < v_near.t
                        v_near = tail
                        # @assert tail.i ≠ head.i # otherwise, we have tail == v_near already
                    end

                    if Ef ≥ Eb
                        λ₋ += ΔE
                    else
                        λ₊ += ΔE
                    end
                    Δt /= (x.β * λ₊)
                    if Δt < t_eps
                        @goto CYCLE_START🔁
                    end
                    t_new = t + Δt
                    t_bound = v_near.t - t_eps
                    if (t_new < t_bound) && metro(δ == 0 ? (λ₋ / λ₊) : inv(λ₊))
                        # no interaction encountered
                        head <<= t_new
                        li[head_id] = head
                        δ = i8(0)
                        loc = _at_free
                        cycle_size += 1
                        # elseif (t_new ≥ t_bound) && metro(δ == 0 ? λ₋ : 1.)
                    elseif (t_new ≥ t_bound) && (δ ≠ 0 || metro(λ₋))
                        # interaction encountered
                        head <<= prevfloat(v_near.t)
                        li[head_id] = head
                        δ = i8(-1)
                        # @assert v_near.j ≠ head.i
                        loc = if v_near.op == I_
                            _at_dummy
                        elseif v_near.j == IndexType(0) # v_near is the tail
                            v_near.i == head.i ? _at_stop : _at_green
                        elseif v_near.i == head.i
                            _at_kink
                        else
                            # @assert v_near.i ≠ head.i && v_near.j ≠ head.i
                            $(
                                if is_KCM_Kagome
                                    quote
                                        kcm_nb_kink = v_near
                                    end
                                end
                            )
                            _at_nbkink
                        end::WormLocation
                        cycle_size += 1
                    end
                else # if D == -1 # move backward <--
                    if v_near.t == 1.0 # dummy element is at both 0 and β, here we regard it at 0.
                        v_near <<= 0.0
                    end
                    @nexprs $znbs k -> begin
                        if assoc_k.i > IndexType(0) && (v_near.t < assoc_k.t < 1.0)
                            v_near = assoc_k
                        end
                    end
                    # however, the head cannot pass the tail
                    if v_near.t < tail.t < head.t
                        v_near = tail
                        # @assert tail.i ≠ head.i # otherwise, we have tail == v_near already
                    end
                    if Ef ≥ Eb
                        λ₋ += ΔE
                    else
                        λ₊ += ΔE
                    end
                    Δt /= (x.β * λ₊)
                    if Δt < t_eps
                        @goto CYCLE_START🔁
                    end
                    t_new = t - Δt
                    t_bound = v_near.t + t_eps
                    if (t_new > t_bound) && metro(δ == 0 ? (λ₋ / λ₊) : inv(λ₊))
                        # no interaction
                        head <<= t_new
                        li[head_id] = head
                        δ = i8(0)
                        loc = _at_free
                        cycle_size += 1
                        # elseif (t_new ≤ t_bound) && metro(δ == 0 ? λ₋ : 1.)
                    elseif (t_new ≤ t_bound) && (δ ≠ 0 || metro(λ₋))
                        head <<= nextfloat(v_near.t)
                        li[head_id] = head
                        δ = i8(+1)
                        loc = if v_near.op == I_
                            _at_dummy
                        elseif v_near.j == IndexType(0) # v_near is the tail
                            v_near.i == head.i ? _at_stop : _at_green
                        elseif v_near.i == head.i
                            _at_kink
                        else
                            # @assert v_near.i ≠ head.i && v_near.j ≠ head.i
                            $(
                                if is_KCM_Kagome
                                    quote
                                        kcm_nb_kink = v_near
                                    end
                                end
                            )
                            _at_nbkink
                        end::WormLocation
                        cycle_size += 1
                    end
                end

                @goto CYCLE_START🔁

            elseif dice < F_ins  # [INSERT_KINK]
                # println("   inserting...")
                if loc ≠ _at_free
                    @goto CYCLE_START🔁
                end
                D = randsign()
                j = rand(hps) |> IndexType
                if j == IndexType(0)
                    @goto CYCLE_START🔁
                end
                lj = x[j]
                B = StateType(head.op)
                jqL, jqR = vindex_around(lj, head.t)
                nj = lj[jqR].n_L # @assert nj == lj[jqL].n_R
                nmid = nj - B * D
                if nmid > H.nmax || nmid < i8(0)
                    @goto CYCLE_START🔁
                end
                
                # case with the kinetic constraint, check the third one
                Wk = max(nj, nmid) * bond_weight(H, i, j)
                $(
                    if is_KCM_Kagome
                        quote
                            if slice_at(x[control_site(H, i, j)], head.t) == StateType(0)
                                Wk *= H.ξ
                            end
                        end
                    end
                )
                P_acc = $(2*zhps) * Wk * P_del2ins
                # P_acc = Wk/(Wk + C_kinks)

                if metro(P_acc)
                    hps_j = get_hps(H, j)::NTuple{$zhps,Int}
                    @nexprs $zhps k -> begin
                        if hps_j[k] > IndexType(0) && hps_j[k] ≠ i && close_to_any(x[hps_j[k]], head.t) # then insertion fail
                            @goto CYCLE_START🔁
                        end
                    end
                    
                    hps = hps_j
                    Kj = Element(head.t, j, i,
                        D == i8(+1) ? nj : nmid,
                        D == i8(+1) ? nmid : nj,
                        OpType(-B)
                    )
                    li[head_id] = Element(head.t, head.i, j, head.n_L, head.n_R, head.op)
                    head = Element(nextfloat(head.t, D), j, IndexType(0), Kj.n_R, Kj.n_L, head.op)
                    # update worm
                    if D == i8(+1)
                        insert!(lj, jqR, Kj)
                        insert!(lj, jqR + 1, head)
                        head_id = jqR + 1
                    else
                        insert!(lj, jqR, head)
                        insert!(lj, jqR + 1, Kj)
                        head_id = jqR
                    end
                    δ = D
                    loc = _at_kink
                    cycle_size += 1
                    i = head.i
                    nbs = get_nbs(H, head.i)
                    li = lj
                end
                @goto CYCLE_START🔁

            elseif dice < F_del  # [DELETE_KINK]
                # println("deleting...")
                if loc ≠ _at_kink
                    @goto CYCLE_START🔁
                end
                Ki_id = head_id - δ
                Ki = li[Ki_id]
                # @assert head.t == nextfloat(Ki.t, δ)
                if head.op == Ki.op # case 1 : pass kink with prob 1
                    @assert H.nmax > StateType(1)
                    li[head_id] = Element(Ki.t,
                        Ki.i, Ki.j, head.n_L, head.n_R, Ki.op
                    ) # now head_id is Ki
                    li[Ki_id] = head = Element(nextfloat(Ki.t, -δ),
                        head.i, head.j, Ki.n_L, Ki.n_R, head.op
                    )
                    i = head.i
                    li = li
                    head_id = Ki_id
                    δ = -δ
                    loc = _at_kink
                    cycle_size += 1
                else # case 2 : delete_kink
                    j = Ki.j
                    lj = x[j]
                    Kj_id = vindex(lj, Ki.t)
                    Kj = lj[Kj_id]
                    # @assert Kj.t == Ki.t
                    Wk = max(head.n_L, head.n_R) * bond_weight(H, i, j)
                    $(
                        if is_KCM_Kagome
                            quote
                                if slice_at(x[control_site(H, i, j)], head.t) == StateType(0)
                                    Wk *= H.ξ
                                end
                            end
                        end
                    )
                    P_acc = inv($(2*zhps) * Wk * P_del2ins)
                    # P_acc = C_kinks / (Wk + C_kinks)
                    if metro(P_acc)
                        head = Element(Kj.t, Kj.i, IndexType(0), Kj.n_L, Kj.n_R, head.op)
                        lj[Kj_id] = head
                        del_id = min(Ki_id, head_id)
                        deleteat!(li, range(del_id, length=2))
                        # update worm
                        i = head.i
                        li = lj
                        head_id = Kj_id
                        hps = get_hps(H, head.i)
                        nbs = get_nbs(H, head.i)
                        δ = i8(0)
                        loc = _at_free
                        cycle_size += 1
                    end
                end
                @goto CYCLE_START🔁

            elseif dice < F_glue # [GLUE_WORM]
                if loc == _at_stop && metro(inv(2 * Cw * max(head.n_L, head.n_R) * P_glue))
                    # @assert head.i == tail.i && head.t == nextfloat(tail.t, δ)
                    del_id = min(head_id, head_id - δ)
                    deleteat!(li, range(del_id, length=2))
                    # @assert check_wl(li)
                    cycle_size += 1
                    continue
                else
                    @goto CYCLE_START🔁
                end
            else # update boson
                $(
                    if is_Bond_Boson
                        quote
                            update_rand_boson!(H, x, fB)
                            cycle_size += 1
                            @goto CYCLE_START🔁
                        end
                    end
                )
            end
        end
        return cycle_size
    end
end