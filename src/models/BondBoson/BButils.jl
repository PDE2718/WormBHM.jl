using Random

struct IntRangeDist
    _bias::Int
    cdf::Vector{Float64}
end
function (f::IntRangeDist)(rng=Random.default_rng())::Int
    return searchsortedfirst(f.cdf, rand(rng)) + f._bias
end
function lims(f::IntRangeDist)::NTuple{2,Int}
    return f._bias .+ (firstindex(f.cdf), lastindex(f.cdf))
end

function _ratio_powk(k, np, n0)
    @assert k ≥ 0 && np ≥ 0 && n0 ≥ 0
    return if k == 0
        1.0
    elseif np == 0 || n0 == 0
        0.0
    else
        (np / n0)^k
    end
end

function bb_dist_k(k::Int, a::Float64, b::Float64, nBmax::Int)
    @assert k ≥ 0 && nBmax > 0
    if (a > 0.0 || (a == 0.0 && b > 0.0))

    else
        println("a = $(a), b=$(b)")
        error()
    end
    P = Float64[1.0]
    nl = round(
        a > 0.0 ? ((-b + √(b^2 + 8a * k)) / (4a)) : (k / b),
        RoundDown
    )
    nl = clamp(nl, 1, nBmax)
    nr = nl
    while P[end] > eps(Float64) && nr < nBmax
        n0 = nr
        np = nr + 1
        Pnext = P[end] * _ratio_powk(k, np, n0) * exp((a * (np + n0) + b) * (n0 - np))
        push!(P, Pnext)
        nr = np
    end
    while P[1] > eps(Float64) && nl > 0
        n0 = nl
        np = nl - 1
        Pprev = P[1] * _ratio_powk(k, np, n0) * exp((a * (np + n0) + b) * (n0 - np))
        pushfirst!(P, Pprev)
        nl = np
    end
    P = cumsum(P)
    P ./= P[end]
    while P[end] > prevfloat(1.0)
        pop!(P)
    end
    while P[1] < eps(Float64)
        popfirst!(P)
        nl += 1
    end
    dist = IntRangeDist(nl - 1, P)
    dist_lims = lims(dist)
    # @assert 0 ≤ dist_lims[1] < dist_lims[2] ≤ nBmax
    if !(0 ≤ dist_lims[1] < dist_lims[2] ≤ nBmax)
        println(dist)
        for i ∈ 1:10
            println(dist())
        end
        println("k = $(k)")
        println(dist_lims)
    end
    return dist
end

# P(n;k) = n^k * exp(- a n² - b n)
struct BondSampler
    a::Float64
    b::Float64
    nBmax::Int
    dists::Vector{IntRangeDist}
end

function BondSampler(β::Float64, μ::Float64, U::Float64, nBmax::Int)
    a = β * 0.5U
    b = β * (-μ - 0.5U)
    dists = IntRangeDist[]
    return BondSampler(a, b, nBmax, dists)
end

function (f::BondSampler)(k::Int, rng=Random.default_rng())::Int
    @assert k ≥ 0
    D = f.dists
    k0 = lastindex(D)
    for ki ∈ k0:k
        push!(D, bb_dist_k(ki, f.a, f.b, f.nBmax))
    end
    sp = D[k+1](rng)
    @assert 0 ≤ sp ≤ f.nBmax
    return sp
end


@inline function _boson_energy_diff(np::Int, n::Int, μ::f64, U::f64)::f64
    return 0.5 * U * (np * (np - 1) - n * (n - 1)) - μ * (np - n)
end

function update_rand_boson!(H::Ham, x::Wsheet{Nw}, fB::Nothing) where {Ham<:BH_Parameters, Nw}
    B = H.bosons
    β = x.β
    B_inds = eachindex(B)
    k = rand(B_inds)
    i, j = bond_sites(H, k)
    M = count(e -> e.j == j, x[i])
    n = B[k]
    np = n + round(Int, randn() / (2β), RoundFromZero)
    if 0 ≤ np ≤ H.nBmax
        Δ = _boson_energy_diff(np, n, H.μb, H.Ub)
        P_acc = ((np^M) / (n^M)) * exp(-β * Δ)
        if metro(P_acc)
            B[k] = np
        end
    end

    kp = k
    while kp == k
        k = rand(B_inds)
        kp = rand(B_inds)
    end

    i, j = bond_sites(H, k)
    M = count(e -> e.j == j, x[i])
    n = B[k]

    i, j = bond_sites(H, kp)
    Mp = count(e -> e.j == j, x[i])
    np = B[kp]

    if np ≠ n
        W1 = n^M * np^Mp
        W2 = np^M * n^Mp
        r = (n / np)^(Mp - M)
        if metro(r)
            B[k] = np
            B[kp] = n
        end
    end
    return nothing
end

function update_rand_boson!(H::Ham, x::Wsheet{Nw}, fB::BondSampler) where {Ham<:BH_Parameters,Nw}
    B = H.bosons
    B_inds = eachindex(B)
    k = rand(B_inds)
    i, j = bond_sites(H, k)
    M = count(e -> e.j == j, x[i])
    n = fB(M)
    B[k] = n
    kp = rand(B_inds)

    if kp ≠ k
        i, j = bond_sites(H, kp)
        Mp = count(e -> e.j == j, x[i])
        np = B[kp]

        if np ≠ n && n≠0 && np≠0 && metro((n / np)^(Mp - M))
            B[k] = np
            B[kp] = n
        end
    end

    return nothing
end

# function _count_j(l::Wline, j)::Int
#     ESIZE = sizeof(Element)
#     nb = 0
#     tar = i32(j)
#     l0 = (l |> pointer |> Ptr{i32}) + sizeof(f64) + sizeof(i32)
#     @inbounds for p ∈ range(
#         start=l0,
#         step=sizeof(Element),
#         length=length(l),
#     )
#         nb += (unsafe_load(p) == tar)
#     end
#     return nb
# end