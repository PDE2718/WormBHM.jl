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
include("simulation.jl")
include("extra_tools.jl")
include("custom/BBCubic_simulation.jl")


for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end