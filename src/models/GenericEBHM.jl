struct SiteInfo{Nk,Nv}
    μ::f64
    U::f64
    jk::NTuple{Nk,Int}  # changed from hops to jk
    K::NTuple{Nk,f64}   # changed from J to K
    jv::NTuple{Nv,Int}  # changed from nbs to jv
    V::NTuple{Nv,f64}
end
function check_EBHM(para::Array{SiteInfo{Nk,Nv},Nw}) where {Nw,Nk,Nv}
    for (i, para_i) ∈ enumerate(para)
        for (j, val) ∈ zip(para_i.jk, para_i.K)  # changed from hops to jk and J to K
            if j ≠ 0
                @assert count(==(j), para_i.jk) == 1  # changed from hops to jk
                para_j = para[j]
                @assert para_j.K[findfirst(==(i), para_j.jk)] == val  # changed from J to K
            else
                @assert val == 0.0
            end
        end
        for (j, val) ∈ zip(para_i.jv, para_i.V)  # changed from nbs to jv
            if j ≠ 0
                @assert count(==(j), para_i.jv) == 1  # changed from nbs to jv
                para_j = para[j]
                @assert para_j.V[findfirst(==(i), para_j.jv)] == val  # changed from nbs to jv
            else
                @assert val == 0.0
            end
        end
    end
    return true
end
@kwdef struct GenericEBHM{IsLattice,Nw,Nk,Nv} <: BH_Parameters  
    nmax::StateType = 1 # nmax must be implemented
    para::Array{SiteInfo{Nk,Nv},Nw}  
end
function GenericEBHM_lattice(nmax::StateType, para::Array{SiteInfo{Nk,Nv},Nw}) where {Nw,Nk,Nv}
    @assert check_EBHM(para)
    @assert nmax ≥ 1
    return GenericEBHM{true,Nw,Nk,Nv}(nmax, para)
end
function get_hps(H::GenericEBHM, i)
    return H.para[i].jk
end
function get_nbs(H::GenericEBHM, i)
    return H.para[i].jv
end
function N_wldim(::Type{GenericEBHM{IsLattice,Nw,Nk,Nv}}) where {IsLattice,Nw,Nk,Nv}
    return Nw
end
function N_hops(::Type{GenericEBHM{IsLattice,Nw,Nk,Nv}}) where {IsLattice,Nw,Nk,Nv}
    return Nk
end
function N_nbs(::Type{GenericEBHM{IsLattice,Nw,Nk,Nv}}) where {IsLattice,Nw,Nk,Nv}
    return Nv
end
function bond_weight(H::GenericEBHM, i::Integer, j::Integer)::f64
    para_i = H.para[i]
    k::Union{Nothing,Int} = findfirst(==(j), para_i.jk)  # changed from hops to jk
    return isnothing(k) ? 0.0 : para_i.K[k]  # changed from J to K
end
# function sublatt_dim(H::GenericEBHM{true,Nw,Nk,Nv}) where {Nw,Nk,Nv}
#     return size(H.para)[Nw]
# end
function diagE(H::GenericEBHM{IsLattice,Nw,Nk,Nv},
    i::Ti, ni::StateType, njs::NTuple{Nv, StateType})::f64 where {
    IsLattice,Nw,Nk,Nv,Ti
}
    P = H.para[i]
    return ni * (0.5 * P.U * (ni-StateType(1)) + sum(P.V .* njs) - P.μ)
end
# function site_diff(H::GenericEBHM{true,Nw,Nk,Nv}, i::Ti, j::Tj) where {Nw,Nk,Nv,Ti,Tj}
#     L = size(H.para)
#     lattice = L |> CartesianIndices
#     ri = lattice[i] |> Tuple
#     rj = lattice[j] |> Tuple
#     D = mod1.((ri .- rj .+ 1), L)[1:(Nw-1)]
#     return CartesianIndex((D..., ri[Nw], rj[Nw]))
# end

function simple_measure!(
    m::SimpleMeasure,
    x::Wsheet{Nw},
    H::GenericEBHM{IsLattice,Nw,Nk,Nv},
    bond_buffer::Wline
)::Nothing where {IsLattice,Nw,Nk,Nv}
    # display(x)
    μ = U = V = 0.0
    N = 0
    NK = NKp = 0
    lattice = CartesianIndices(x.wl)
    @inbounds for (i, c) ∈ enumerate(lattice)
        li = x[i]
        para_i = H.para[i]
        ni¹, ni² = integrated_density(Val(2), li)
        μ -= (para_i.μ * ni¹)
        U += (0.5 * para_i.U * (ni²-ni¹))
        for (j, Vj) ∈ zip(para_i.jv, para_i.V)
            if 0 < j < i
                V += Vj * integrated_densities(li, x[j])
            end
        end
        for e ∈ li
            if e.op == b_
                NK += 1
            elseif e.op == b̂_
                NKp += 1
            end
        end
    end
    @assert NK == NKp
    K = -NK / x.β
    E = μ + U + V + K
    m.props = m.props .+ (E, abs2(E), N, abs2(N), μ, U, V, K, V_new)
    m.n_measure += 1
    return nothing
end

function simple_measure_names(H::GenericEBHM)
    return (:E, :E², :N, :N², :μ, :U, :V, :K)
end

function Wsheet(β::f64, H::GenericEBHM{IsLattice,Nw,Nk,Nv}) where {IsLattice,Nw,Nk,Nv}
    Wsheet(β, zeros(StateType, size(H.para)))
end
