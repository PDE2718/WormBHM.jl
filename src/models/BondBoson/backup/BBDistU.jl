using Random

struct BBDistU # P(n) = n^k * exp(- a n² - b n)
    k::Int
    a::Float64
    b::Float64
    _bias::Int
    cdf::Vector{Float64}
end
function BBDistU(k::Int, a::Float64, b::Float64)
    @assert k ≥ 0
    @assert (a > 0.0 || (a == 0.0 && b > 0.0))
    P = 1.0
    n = 0
    _bias = 1
    cdf = Float64[1.0]
    S = Snew = 1.0
    while n < 1000000
        n += 1
        P = exp(k * log(n) - a * abs2(n) - b * n)
        Snew += P
        push!(cdf, Snew)
        if S / Snew > prevfloat(1.0)
            break
        else
            S = Snew
        end
    end
    cdf ./= cdf[end]
    @assert issorted(cdf)
    while cdf[end] > prevfloat(1.0)
        pop!(cdf)
    end
    push!(cdf, 1.0)
    while cdf[1] < eps(Float64)
        popfirst!(cdf)
        _bias += 1
    end
    resize!(cdf, length(cdf))
    return BBDistU(k, a, b, _bias, cdf)
end
function (f::BBDistU)(rng=Random.default_rng())
    return searchsortedfirst(f.cdf, rand(rng)) - f._bias
end

function rand_boson!(fB::Vector{BBDist}, Nk::Int, βω::Float64)
    @assert Nk ≥ 0
    @assert βω > 0.0
    l0 = length(fB)
    for k ∈ l0:Nk
        push!(fB, BBDist(k, βω))
    end
    f = fB[Nk+1]
    if f.b == βω && f.k == Nk
        return f()
    else
        empty!(fB)
        return rand_boson!(fB, Nk, βω)
    end
end

function rand_boson!(fB::Vector{BBDistU}, Nk::Int, a::Float64, b::Float64)
    @assert Nk ≥ 0
    l0 = length(fB)
    for k ∈ l0:Nk
        push!(fB, BBDistU(k, a, b))
    end
    f = fB[Nk+1]
    if f.a == a && f.b == b && f.k == Nk
        return f()
    else
        empty!(fB)
        return rand_boson!(fB, Nk, a, b)
    end
end

function update_rand_boson!(fB::Vector{BBDist}, H, x)
    B = H.bosons
    kij = rand(eachindex(B))
    i, j = bond_sites(H, kij)
    Nk = count(e -> e.j == j, x[i])
    np = rand_boson!(fB, Nk, x.β * H.ω)
    if 0 ≤ np ≤ H.nBmax
        B[kij] = np
    end
    return nothing
end

function update_rand_boson!(fB::Vector{BBDistU}, H, x)
    B = H.bosons
    a = x.β * 0.5 * H.Ub
    b = -x.β * (0.5 * H.Ub + H.μb)
    kij = rand(eachindex(B))
    i, j = bond_sites(H, kij)
    Nk = count(e -> e.j == j, x[i])
    np = rand_boson!(fB, Nk, a, b)
    while np > H.nBmax
        np = rand_boson!(fB, Nk, a, b)
    end
    if 0 ≤ np ≤ H.nBmax
        B[kij] = np
    end
    return nothing
end

function update_rand_boson_trivial!(H, x)
    B = H.bosons
    β = x.β
    a = β * 0.5 * H.Ub
    b = -β * (0.5 * H.Ub + H.μb)
    B_inds = eachindex(B)
    k = rand(B_inds)
    i, j = bond_sites(H, k)
    M = count(e -> e.j == j, x[i])
    n = B[k]
    np = n + round(Int, randn()/(2β), RoundFromZero)
    if 0 ≤ np ≤ H.nBmax
        P_acc = ((np^M) / (n^M)) * exp(n*(a*n+b)-np*(a*np+b))
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
        r = (n/np)^(Mp-M)
        if metro(r)
            B[k] = np
            B[kp] = n
        end
    end
    return nothing
end