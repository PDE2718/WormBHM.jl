@generated function simulate!(
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
    sweep_size_limit::Tuple{Int,Int}=(1, 1024),
    shuffle_snapshot::f64=0.0,
    ) where {Nw, Ham <: BH_Parameters,
        T_GreenFuncBin <: Union{GreenFuncBin,Nothing},
        T_SimpleMeasure <: Union{SimpleMeasure, Nothing},
        T_DensityMap <: Union{DensityMap, Nothing},
        T_StructureFactor <: Union{StructureFactorND, Nothing},
        T_Snaps <: Union{FixedSizeRecorder, Nothing}
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
        sweep_size = 100
        while time_ns() < t_limit
            total_cycle_size += worm_cycle!(x, H, 
                update_consts, cycle_probs,
                sweep_size, nothing
            )
            N_cycle += sweep_size
        end
        average_size = total_cycle_size / N_cycle
        Nverts = sum(length, x.wl) - length(x.wl)
        sweep_size = ceil(Int, Nverts / average_size)
        sweep_size = clamp(sweep_size, sweep_size_limit...)
        
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
            total_cycle_size += worm_cycle!(x, H,
                update_consts, cycle_probs,
                sweep_size, Gf
            )
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
