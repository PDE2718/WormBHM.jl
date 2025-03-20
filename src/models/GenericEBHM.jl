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
@kwdef struct GenericEBHM{IsLattice,Nw,Nk,Nv} <: BH_Parameters{IsLattice,Nw,Nk,Nv}
    nmax::StateType = 1 # nmax must be implemented
    para::Array{SiteInfo{Nk,Nv},Nw}  
end
function GenericEBHM_lattice(nmax::StateType, para::Array{SiteInfo{Nk,Nv},Nw}) where {Nw,Nk,Nv}
    @assert check_EBHM(para)
    @assert nmax ≥ 1
    return GenericEBHM{true,Nw,Nk,Nv}(nmax, para)
end
function wl_size(H::GenericEBHM)
    return size(H.para)
end
function get_hps(H::GenericEBHM, i)
    return H.para[i].jk
end
function get_nbs(H::GenericEBHM, i)
    return H.para[i].jv
end

function bond_weight(H::GenericEBHM{IsLattice,Nw,Nk,Nv},
    i::Integer, j::Integer)::f64 where {IsLattice,Nw,Nk,Nv}
    para_i::SiteInfo{Nk, Nv} = H.para[i]
    k::Union{Nothing,Int} = findfirst(==(j), para_i.jk)  # changed from hops to jk
    return isnothing(k) ? 0.0 : para_i.K[k]  # changed from J to K
end

function diagE(H::GenericEBHM{IsLattice,Nw,Nk,Nv},
    i::Ti, ni::StateType, njs::NTuple{Nv, StateType})::f64 where {
    IsLattice,Nw,Nk,Nv,Ti
}
    P = H.para[i]
    return ni * (0.5 * P.U * (ni-StateType(1)) + sum(P.V .* njs) - P.μ)
end
function simple_measure_names(H::Type{T}) where {T<:GenericEBHM}
    return (:E, :E², :N, :N², :μ, :U, :V, :K)
end
function simple_measure!(
    m::SimpleMeasure{8},
    x::Wsheet{Nw},
    H::GenericEBHM{IsLattice,Nw,Nk,Nv},
)::Nothing where {IsLattice,Nw,Nk,Nv}
    μ = U = V = 0.0
    N = 0.0
    NK = NKp = 0
    lattice = CartesianIndices(x.wl)
    for (i, c) ∈ enumerate(lattice)
        li = x[i]
        N += li[1].n_L
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
    m.props = m.props .+ (E, abs2(E), N, abs2(N), μ, U, V, K)
    m.n_measure += 1
    return nothing
end

function set_Kij!(para::Array{SiteInfo{Nk,Nv},Nw}, 
    i::Int, j::Int, K) where {Nw,Nk,Nv}

    @assert checkbounds(Bool, para, i) && checkbounds(Bool, para, j)
    site_info_size::Int = 2 + 2Nk + 2Nv
    x = reinterpret(Int, para)
    st = site_info_size * (i - 1) + 2
    Kval = reinterpret(Int, K)
    @inbounds for q ∈ (st+1):(st+Nk)
        if x[q] == j
            x[q+Nk] = Kval
            break
        end
    end
    i, j = (j, i)
    st = site_info_size * (i - 1) + 2
    @inbounds for q ∈ (st+1):(st+Nk)
        if x[q] == j
            x[q+Nk] = Kval
            break
        end
    end
    return nothing
end
function set_Vij!(para::Array{SiteInfo{Nk,Nv},Nw},
    i::Int, j::Int, V) where {Nw,Nk,Nv}

    @assert checkbounds(Bool, para, i) && checkbounds(Bool, para, j)
    site_info_size::Int = 2 + 2Nk + 2Nv
    x = reinterpret(Int, para)
    st = site_info_size * (i - 1) + 2 + 2Nk
    Vval = reinterpret(Int, V)
    @inbounds for q ∈ (st+1):(st+Nv)
        if x[q] == j
            x[q+Nv] = Vval
            break
        end
    end
    i, j = (j, i)
    st = site_info_size * (i - 1) + 2 + 2Nk
    @inbounds for q ∈ (st+1):(st+Nv)
        if x[q] == j
            x[q+Nv] = Vval
            break
        end
    end
    return nothing
end
function set_μi!(para::Array{SiteInfo{Nk,Nv},Nw},
    i::Int, μ::Float64) where {Nw,Nk,Nv}

    reinterpret(Float64, para)[(2+2Nk+2Nv)*(i-1)+1] = μ
    return nothing
end
function set_Ui!(para::Array{SiteInfo{Nk,Nv},Nw},
    i::Int, U::Float64) where {Nw,Nk,Nv}

    reinterpret(Float64, para)[(2+2Nk+2Nv)*(i-1)+2] = U
    return nothing
end
function bond_iter_K(para::Array{SiteInfo{Nk,Nv},Nw}) where {Nw,Nk,Nv}
    return ((i, j) => K for (i, p) ∈ enumerate(para) for (j, K) ∈ zip(p.jk, p.K) if (i < j))
end
function bond_iter_V(para::Array{SiteInfo{Nk,Nv},Nw}) where {Nw,Nk,Nv}
    return ((i, j) => V for (i, p) ∈ enumerate(para) for (j, V) ∈ zip(p.jv, p.V) if (i < j))
end