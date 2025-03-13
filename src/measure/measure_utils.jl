@generated function integrated_density(::Val{N}, l::Wline)::NTuple{N, f64} where {N}
    quote
        # @assert issorted(l)
        t::f64 = 0.0
        n::StateType = l[1].n_L
        @nexprs $(N) i->begin
            n_i::f64 = 0.
        end
        @inbounds for e ∈ l
            Δt = e.t - t
            np::Int = 1
            @nexprs $(N) i -> begin
                np *= e.n_L
                n_i += (Δt * np)
            end
            t = e.t
        end
        res = @ntuple $(N) i->n_i
        return res
    end
end

@generated function integrated_densities(vs::Vararg{Wline,N})::f64 where {N}
    @assert N > 1 "for one world line, use `integrated_density` instead"
    quote
        L0 = 0
        @nexprs $(N) i -> begin
            v_i = vs[i]
            L_i = length(v_i)
            t_i = 1
            n_i = v_i[1].n_L
            L0 += L_i
        end
        t = ∏ = 0.0
        @inbounds for t_dst ∈ 1:L0
            min_idx = 0
            local min_val::Element
            @nexprs $(N) i -> begin
                if t_i ≤ L_i && (min_idx == 0 || v_i[t_i] < min_val)
                    min_val = v_i[t_i]
                    min_idx = i
                end
            end
            if min_val.t > t
                ∏ += (min_val.t - t) * prod(
                    @ntuple $(N) i -> n_i
                )
                t = min_val.t
            end
            @nexprs $(N) i -> if i == min_idx
                t_i += 1
                n_i = min_val.n_R
            end
        end
        return ∏
    end
end

@generated function merge_sorted!(dst::Vector{T}, vs::Vararg{Vector{T},N})::Nothing where {T,N}
    if N == 1
        quote
            v = first(vs)
            resize!(dst, length(v))
            copy!(dst, v)
            return nothing
        end
    else
        quote
            L0 = 0
            @nexprs $(N) i -> begin
                v_i = vs[i]
                L_i = length(v_i)
                t_i = 1
                L0 += L_i
            end
            resize!(dst, L0)
            @inbounds for t_dst ∈ 1:L0
                min_idx = 0
                local min_val::T
                @nexprs $(N) i -> begin
                    if t_i ≤ L_i && (min_idx == 0 || v_i[t_i] < min_val)
                        min_val = v_i[t_i]
                        min_idx = i
                    end
                end
                dst[t_dst] = min_val
                @nexprs $(N) i -> if i == min_idx
                    t_i += 1
                end
            end
            return nothing
        end
    end
end