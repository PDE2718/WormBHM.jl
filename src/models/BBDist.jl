using Random
struct BBDist # P(n) = n^k * exp(-b*n)
    k::Int
    b::Float64
    _bias::Int
    cdf::Vector{Float64}
end
function BBDist(k::Int, b::Float64)
    @assert k ≥ 0
    @assert b > 0.0
    P = 1.0
    n = 0
    _bias = 1
    cdf = Float64[1.0]
    S = Snew = 1.0
    while n < 1000000
        n += 1
        P = exp(k * log(n) - b * n)
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
    return BBDist(k, b, _bias, cdf)
end

function (f::BBDist)(rng=Random.default_rng())
    return searchsortedfirst(f.cdf, rand(rng)) - f._bias
end

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
    if f.b == βω
        return fB[Nk+1]()
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
    if f.a == a && f.b == b
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
    if 0 ≤ np ≤ H.nBmax
        B[kij] = np
    end
    return nothing
end