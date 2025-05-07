@inline function _boson_energy_diff(np::Int, n::Int, μ::f64, U::f64)::f64
    return 0.5 * U * (np * (np - 1) - n * (n - 1)) - μ * (np - n)
end

function update_rand_boson!(H::Ham, x::Wsheet{Nw}) where {Ham<:BH_Parameters, Nw}
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