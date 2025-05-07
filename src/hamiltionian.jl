abstract type BH_Parameters{IsLattice, Nw, Nk, Nv} end

islattice(::Type{<:BH_Parameters{IsLattice,Nw,Nk,Nv}}) where {
    IsLattice,Nw,Nk,Nv
} = IsLattice
N_wldim(::Type{<:BH_Parameters{IsLattice,Nw,Nk,Nv}}) where {
    IsLattice,Nw,Nk,Nv
} = Nw
N_hps(::Type{<:BH_Parameters{IsLattice,Nw,Nk,Nv}}) where {
    IsLattice,Nw,Nk,Nv
} = Nk
N_nbs(::Type{<:BH_Parameters{IsLattice,Nw,Nk,Nv}}) where {
    IsLattice,Nw,Nk,Nv
} = Nv

function wl_size end
function get_hps end
function get_nbs end

function diagE end
function bond_weight end
function simple_measure_names end
function simple_measure! end

function Wsheet(β::f64, H::T) where {T<:BH_Parameters}
    @assert β > 0
    Wsheet(β,[empty_wline(i, StateType(0)) for i ∈ LinearIndices(wl_size(H))])
end

function simple_measure!(::Nothing, x, H)
    return nothing
end
function simple_measure_count(H::Type{T}) where {T<:BH_Parameters}
    return simple_measure_names(H) |> length
end
function simple_measure_names(::T) where {T<:BH_Parameters}
    return simple_measure_names(T)
end
function simple_measure_count(::T) where {T<:BH_Parameters}
    return simple_measure_count(T)
end
function SimpleMeasure(H::Type{Ham}) where {Ham<:BH_Parameters}
    return SimpleMeasure(
        ntuple(x -> 0.0, simple_measure_count(H)),
        simple_measure_names(H),
        0
    )
end
function SimpleMeasure(::T) where {T<:BH_Parameters}
    return SimpleMeasure(T)
end
function init_wsheet(β::T_β, H::T_H) where {T_β<:Real,T_H<:BH_Parameters}
    @assert β > T_β(0)
    return Wsheet(Float64(β), zeros(StateType, wl_size(H)))
end
include("models/BH_Square.jl")
include("models/BH_Cubic.jl")
include("models/Kagome.jl")
include("models/GenericEBHM.jl")
include("models/HybridKagome.jl")
include("models/BondBoson/BButils.jl")
include("models/BondBoson/BBChainU.jl")
include("models/BondBoson/BBSquareU.jl")
include("models/BondBoson/BBCubicU.jl")
##########################################################################################
