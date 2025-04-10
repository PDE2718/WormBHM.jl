@generated function simulate_bb!(
    x::Wsheet{Nw},
    H::Ham,
    update_consts::UpdateConsts,
    cycle_probs::CycleAccumProb,
    time_ther::Second,
    time_simu::Second,
    simple_measure::T_SimpleMeasure=nothing,
    density_map::T_DensityMap=nothing,
    Gf::T_GreenFuncBin=nothing,
    Sk::T_StructureFactor=nothing,
    snapshots::T_Snaps=nothing,
    bosonhist::T_BosonHist=nothing,
    density_map_boson::T_DensityMapBoson=nothing,
    optimal_update_consts=true,
    optimal_cycle_probs=true,
    sweep_size_limit::Tuple{Int,Int}=(1, 1024),
    shuffle_snapshot::f64=0.0,
    ) where {Nw, Ham <: BH_Parameters,
        T_GreenFuncBin <: Union{GreenFuncBin,Nothing},
        T_SimpleMeasure <: Union{SimpleMeasure, Nothing},
        T_DensityMap <: Union{DensityMap, Nothing},
        T_StructureFactor <: Union{StructureFactorND, Nothing},
        T_Snaps <: Union{FixedSizeRecorder, Nothing},
        T_BosonHist <: Union{Vector{Int}, Nothing},
        T_DensityMapBoson <: Union{DensityMap, Nothing},
    }
    if T_DensityMap == Nothing
        @assert T_StructureFactor == T_Snaps == Nothing "StructureFactor and Snapshots have DensityMap as dependency"
    end
    return quote
        $(
            if T_GreenFuncBin ≠ Nothing
                quote
                    @assert Gf.Cw == update_consts.Cw
                end
            end
        )
        $(
            if T_StructureFactor ≠ Nothing
                quote
                    Sk_plan = rfftplan(Sk)
                end
            end
        )
        @info "thermalizing for $(time_ther)"
        @assert iszero(time_ns() >> 62) "current time_ns is close to wrap and can produce bugs!"
        T_ther::UInt64 = Nanosecond(time_ther).value
        T_simu::UInt64 = Nanosecond(time_simu).value

        t_limit::UInt64 = time_ns() + T_ther
        
        N_cycle = total_cycle_size = 0
        sweep_size = 10
        fB = BBDistBuffer(H)
        while time_ns() < t_limit
            total_cycle_size += worm_cycle_bb!(x, H, 
                update_consts, cycle_probs,
                sweep_size, nothing, fB
            )
            # update_bosons!(H, x, fB)
            N_cycle += sweep_size
        end
        average_size = total_cycle_size / N_cycle
        Nverts = sum(length, x.wl) - length(x.wl)
        sweep_size = ceil(Int, Nverts / average_size)
        sweep_size = clamp(sweep_size, sweep_size_limit...)

        if optimal_update_consts
            @reset update_consts.Eoff = (Nverts/length(x.wl) + 1)/2x.β
            @assert update_consts.Eoff > 0
            @info "reset update_consts.Eoff to optimal : $(update_consts)"
        end
        if optimal_cycle_probs
            Wkbar::Float64 = 1.0 + (H.J * sum(H.bosons) / length(H.bosons))
            # P_acc = $(2*zhps) * Wk * update_consts.P_del2ins
            P_del2ins_opt = inv(Wkbar)
            P_ins = 0.50 * inv(1+P_del2ins_opt)
            P_del = 0.5 - P_ins
            @reset update_consts.P_del2ins = P_del / P_ins
            cycle_probs = CycleProb(0.25, P_ins, P_del, 0.25)
            @info "reset cycle_probs to optimal : $(cycle_probs)"
            @info "reset update_consts to optimal : $(update_consts)"
        end

        # sweep_size = 10
        @assert all(issorted, x.wl)
        @info """Thermalization Statistics
            [total_cycle_size] = $(total_cycle_size)
            [N_cycles        ] = $(N_cycle)
            [average_size    ] = $(average_size)
            [wl length / β   ] = $((mean(length, x.wl)-1) / x.β)
            [N cycle per mes ] = $(sweep_size)
            Now simulating for $(time_simu)
            """

        t_limit = time_ns() + T_simu
        N_cycle = total_cycle_size = 0
        T_measure = N_measure = 0
        while time_ns() < t_limit
            total_cycle_size += worm_cycle_bb!(x, H,
                update_consts, cycle_probs,
                sweep_size, Gf, fB
            )
            # update_bosons!(H,x,fB)
            N_cycle += sweep_size
            tic = time_ns()
            $(
                if T_SimpleMeasure ≠ Nothing
                    quote
                        simple_measure!(simple_measure, x, H)
                    end
                end
            )
            $(
                if T_DensityMap ≠ Nothing
                    quote
                        snapshot!(density_map, x, shuffle_snapshot)
                    end
                end
            )
            $(
                if T_StructureFactor ≠ Nothing
                    quote
                        measure_Sk!(Sk, density_map.ψ, Sk_plan)
                    end
                end
            )
            $(
                if T_Snaps ≠ Nothing
                    quote
                        record!(snapshots, density_map.ψ)
                    end
                end
            )
            $(
                if T_BosonHist ≠ Nothing
                    quote
                        for nb ∈ H.bosons
                            r = nb + 1
                            if length(bosonhist) < r
                                l0 = length(bosonhist)
                                resize!(bosonhist, r)
                                bosonhist[(l0+1):r] .= 0
                            end
                            bosonhist[r] += 1
                        end
                    end
                end
            )
            $(
                if T_DensityMapBoson ≠ Nothing
                    quote
                        density_map_boson.ρ .+= H.bosons
                        density_map_boson.n_measure += 1
                    end
                end
            )
            T_measure += (time_ns() - tic)
            N_measure += 1
        end
        @assert all(issorted, x.wl)
        average_size = total_cycle_size / N_cycle
        @info """[FINISHED] Simulation Statistics
        [total_cycle_size ] = $(total_cycle_size)
        [N_cycles         ] = $(N_cycle)
        [average_size     ] = $(average_size)
        [wl length / β    ] = $((mean(length, x.wl)-1) / x.β)
        [total measure    ] = $(N_measure)
        [total   time (ns)] = $(T_simu)
        [measure time (ns)] = $(T_measure)
        [T_mes / T_tot    ] = $(T_measure/T_simu)
        [SimpleMeasure] : $(simple_measure)
        """
        return nothing
    end
end
