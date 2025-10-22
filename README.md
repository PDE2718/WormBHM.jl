# WormBHM.jl

A Julia implementation of the worm algorithm quantum Monte-Carlo (QMC). Largely inspired by [the original program by Nicolas Sadoune & Lode Pollet](https://github.com/LodePollet/worm).

## Overview
> Warning⚠️: The code is still under development.

`WormBHM.jl` can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:
```bash
pkg> add https://github.com/PDE2718/WormBHM.jl
```

## Minimum example

Consider the example of a soft-core extended Bose-Hubbard model on 2D square lattice, whose hamiltonian reads:
$$
    H = \left(\sum_{\langle i j \rangle} - J (a_i a_j^\dagger + \text{h.c.}) + Vn_in_j \right) + \left( \sum_i  \frac{U}{2} n_i (n_i - 1) - \mu n_i \right)
$$

The parameters set
1. system size $L_x=L_y=8$.
1. inverse temperature $\beta = 8$.
2. periodic boundary condition (PBC).
3. occupation cut-off: $n_{\max}=2$.
4. $J=1.0, V = 0.25, \mu = 1.1, U = 0.5$

We begin by including necessary packages:
```julia
using WormBHM, Dates, Statistics, Logging
```
Now we can define the model hamiltonian and temperature:
```julia
H = BH_Square(nmax=2, Lx=4, Ly=4, U=0.5, J=1.0, V = 0.25, μ=1.1)
β = 8.0 # β = 1/T
hyperpara = UpdateHyperPara(
    Cw=2.0,
    Eoff=1.0,
    P_move=1.0,
    P_ins=1.0,
    P_del=1.0,
    P_glue=1.0,
) # This hyperparameter can be fine tuned to optimize the efficiency.
```
Next, we initialize the world-line configuration:
```julia
x = Wsheet(β, H) # The empty world-line.
```
and we specify the measurements:
```julia
m = (
    simple_measure=SimpleMeasure(H), # basic quantities
    density_map=DensityMap(x),       # particle density histogram 
    Gf=GreenFuncBin(x, hyperpara;    # Green's function
        IsLattice=true,              # whether apply lattice symmetry
        FullImaginary=nothing,       # don't record full imaginary time info.
        lmax=0,                      # set to 0 if the above is `nothing`
        ),
    Sk=StructureFactorND(x),         # Density Structure Factor
    snapshots=SnapShots(x, 128),     # record no more than 128 snapshots
)
simu_opts = (
    optimize_hyperpara=false,   # disable the 
    sweep_size_limit=(1, 1024), # limit of the number of sweeps between measurements.
    shuffle_snapshot=0.0, # [0.0,1.0) for fixed position, 1.0 for random position.
)
```






```julia
using WormQMC
using LinearAlgebra, Statistics, Dates, DelimitedFiles, Logging
# Logging.disable_logging(Logging.Info)

# Example 8×8 hard-core Hubbard
# U ≫ 1 and nmax ≫ 1 can also be set. Here U is useless.
H = BH_Square(nmax=1, Lx=8, Ly=8, U=40.0, J=1.0, V=0.25, μ=0.0)
β = 8.0 # T = 1/β

# Update constants. Can be fine tuned.
update_const = UpdateConsts(0.5, 2.0, 1.0)
cycle_prob = CycleProb(1, 1, 1, 1)

# Thermalization and simulation time. All need to be in second.
time_ther = 2 |> Minute |> Second
time_simu = 20 |> Minute |> Second

# initialize the world line config and measurement
# green_lmax for imaginary-time green's function precision
x = Wsheet(β, H)
m = WormMeasure(x, update_const; green_lmax=100)

#! Do the simulation (should finish in time_ther+time_simu)
onesimu!(x, H, m, update_const, cycle_prob, time_ther, time_simu)

# calculate the density matrix with Gfunc buffer
G0 = normalize_density_matrix(m.Gfunc)[:, :, 1, 1]
# calculate the real-space correlation Cij = ⟨ninj⟩
Cij = Cab(m.Sfact, 1, 1)

# Now display the result
begin
    @info "QMC result"
    m.simple |> display
    m.winding |> display
    println("┌ DensMat")
    writedlm(stdout, G0)
    println("┌ DensCor")
    writedlm(stdout, Cij)
end

# Plot the imaginary green's function
using Plots,LaTeXStrings
τgrid = LinRange(0.0, β, 100)
G0τ = cal_Gτ(m.Gfunc, τgrid)[1]
let G00 = G0τ[1] # normalization with density matrix
    G0τ .*= (G0[1,1] ./ G00)
end
plot(τgrid,G0τ,
    xlims=(0.,β),ylims=(0,0.8),
    xlabel=L"τ",ylabel=L"G(τ,0)",
    framestyle=:box, label = false
)

# check that for hard-core bosons, G(0⁺) + G(0⁻) == 1
@assert ≈(G0τ[1] + G0τ[end], 1, atol=0.05)
```
