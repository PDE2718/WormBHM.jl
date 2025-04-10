module WormBHM
using LinearAlgebra, StaticArrays, Random, Distributions, Statistics, FFTW, LegendrePolynomials, Dates, Accessors
using Base.Cartesian.Base:@ntuple, @nexprs, @nextract

include("element.jl")
include("parameters.jl")
begin
    include("measure/measure_utils.jl")
    include("measure/simple_measure.jl")
    include("measure/density_map.jl")
    include("measure/structure_factor.jl")
    include("measure/greens_function.jl")
    include("measure/snapshots.jl")
end
include("hamiltionian.jl")
include("update.jl")
include("custom/update_bb.jl")
include("simulation.jl")
include("custom/simulation_bb.jl")
include("custom/simulation_opt.jl")
include("extra_tools.jl")

for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end