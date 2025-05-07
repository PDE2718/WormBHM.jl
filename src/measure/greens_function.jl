# Ref: https://ar5iv.labs.arxiv.org/html/1104.3215
# The idea : Legendre Polynomial Expansion
# Gτ = β⁻¹ ∑ₗ P̄ₗ(x) Gₗ
# Gₗ = ∫₀ᵝ dt P̄ₗ(x) Gτ
# P̄(x) = √(2l+1) P(x) : 
#   Legendre Polynomials with schmidt normalization, P(x) is the standard
# here: x(τ) = 2τ/β - 1   ∈ [-1,1]

# accumulator of Matsubara Green Functions & density matrix.
# if not is_full, then measure only the on site Green's function G(τ)
# only for systems with translational symmetry
@inline function xτmap(Δτ::f64, β::f64)::f64
    return (2Δτ / β) - 1.0
end
# collectPl!(G._Pl, x, norm=Val(:schmidt)) # norm with √(2l+1)
function Pl_schmidt_lazy!(_Pl::Vector{Float64}, x::Float64)
    @inbounds begin
        if lastindex(_Pl) ≥ 2 && _Pl[1] == 1.0 && _Pl[2] == x * √3
            return _Pl
        else
            collectPl!(_Pl, x, norm=Val(:schmidt))
            return _Pl
        end
    end
end
function site_diff(IsLattice::Val{true}, x::Wsheet{N}, i::Tid, j::Tid) where {N,Tid<:Integer}
    spatial_shape::NTuple{N - 1,Int}, N_sub::Int = size(x.wl) |> peel_last
    lattice = CartesianIndices(x.wl)
    ri::NTuple{N - 1,Int}, si::Int = lattice[i] |> Tuple |> peel_last
    rj::NTuple{N - 1,Int}, sj::Int = lattice[j] |> Tuple |> peel_last
    dr::NTuple{N - 1,Int} = @. mod1((ri - rj + 1), spatial_shape)
    return CartesianIndex(dr..., si, sj)
end
function site_diff(IsLattice::Val{false}, x::Wsheet, i::Integer, j::Integer)
    @assert checkbounds(Bool, x.wl, i)
    @assert checkbounds(Bool, x.wl, j)
    return CartesianIndex(i, j)
end
mutable struct GreenFuncBin{IsLattice, FullImaginary,
    T_Pl<:Union{Nothing, Vector{f64}},
    T_G0<:Array{Complex{Int}},
    T_Gl<:Union{Nothing, Array{Vector{f64}}}
    }
    const β::f64
    const Cw::f64
    const _Pl::T_Pl # Pₗ buffer for speed
    const G0::T_G0
    const Gl::T_Gl
    insertion_trial::Int
end
function GreenFuncBin(x::Wsheet, hyperpara;
    IsLattice::Bool = true,
    FullImaginary::Union{Nothing, Bool} = nothing,
    lmax::Int = 0,
    )
    @assert x.β > 0.
    Cw = hyperpara.Cw
    @assert Cw > 0.
    wlshape = size(x.wl)
    _Pl = if FullImaginary == nothing
        @assert lmax == 0
        nothing
    else
        @assert lmax > 1
        zeros(lmax + 1)
    end
    G0 = if IsLattice
        @assert length(wlshape) > 1
        spatial_shape = wlshape[1:end-1]
        Nsub = wlshape[end]
        zeros(Complex{Int}, (spatial_shape..., Nsub, Nsub))
    else
        Nsite = length(x.wl)
        zeros(Complex{Int}, (Nsite,Nsite))
    end
    Gl = if FullImaginary == nothing
        nothing
    elseif FullImaginary == false
        [zeros(lmax + 1) for l ∈ x.wl]
    else #if FullImaginary == true
        [zeros(lmax + 1) for g ∈ G0]
    end
    return GreenFuncBin{IsLattice, FullImaginary, typeof(_Pl), typeof(G0), typeof(Gl)}(
        x.β, Cw, _Pl, G0, Gl, 0
    )
end
function Base.show(io::IO,
    G::GreenFuncBin{IsLattice,FullImaginary,T_Pl,T_G0,T_Gl}
    ) where {IsLattice,FullImaginary,T_Pl,T_G0,T_Gl}
    println(io, "┌ GreenFuncBin:")
    println(io, "│ IsLattice = $(IsLattice), FullImaginary = $(FullImaginary)")
    println(io, "│ β = $(G.β), Cw = $(G.Cw)")
    println(io, "│ Pl → $(summary(G._Pl))")
    println(io, "│ G0 → $(summary(G.G0))")
    println(io, "└ Gl → $(summary(G.Gl))")
end

@generated function accumulate!(
    G::GreenFuncBin{IsLattice,FullImaginary,T_Pl,T_G0,T_Gl},
    x::Wsheet{N},
    b::Element, b̂::Element, loc::WormLocation
)::Nothing where {IsLattice,FullImaginary,T_Pl,T_G0,T_Gl,N}
    quote
        @inbounds begin
            if b.op == b̂_
                c = b
                b = b̂
                b̂ = c
            end
            @assert b.op == b_ && b̂.op == b̂_
            Gid = site_diff($(Val(IsLattice)), x, b.i, b̂.i)
            #measure density matrix
            if loc == _at_green
                G.G0[Gid] += complex(1, 0)
            elseif loc == _at_stop
                G.G0[Gid] += (b̂.t > b.t) ? complex(1, 0) : complex(0, 1)
            elseif loc == _at_free 
                $(
                    if FullImaginary == true
                        quote
                            x = 2mod(b̂.t - b.t, 1.0) - 1.0
                            Pl_schmidt_lazy!(G._Pl, x)
                            G.Gl[Gid] .+= G._Pl
                        end
                    elseif FullImaginary == false
                        quote
                            if b̂.i == b.i
                                x = 2mod(b̂.t - b.t, 1.0) - 1.0
                                Pl_schmidt_lazy!(G._Pl, x)
                                G.Gl[b.i] .+= G._Pl
                            end
                        end
                    elseif FullImaginary == nothing
                        @assert T_Gl == Nothing
                        quote
                            nothing
                        end
                    end
                )
            end
            return nothing
        end
    end
end

function cal_Gτ(G::GreenFuncBin, τgrid, lmax::Int=-1)
    @assert G.Gl ≠ nothing "cannot calculate imaginary time GF for bins without Gl"
    Gτs = [zeros(length(τgrid)) for i ∈ CartesianIndices(G.Gl)]
    Plx::Vector{f64} = G._Pl
    lmax = lmax < 0 ? (length(Plx) - 1) : clamp(lmax, 0, length(Plx) - 1)
    for (iτ, τ) ∈ enumerate(τgrid)
        collectPl!(Plx, xτmap(τ, G.β), norm=Val(:schmidt)) # now they are order 0...lmax
        for (Gτ, Gl) ∈ zip(Gτs, G.Gl)
            @inbounds for lp1 ∈ 1:(lmax+1)
                Gτ[iτ] += Gl[lp1] * Plx[lp1]
            end
        end
    end
    Γ = inv(4 * G.Cw * G.insertion_trial)
    for Gτ ∈ Gτs
        Gτ .*= Γ
    end
    return Gτs
end

# two ways of normalizing greens function, as a cross check
function _norm_coeff(G::GreenFuncBin{IsLattice,FullImaginary,T_Pl,T_G0,T_Gl}) where {IsLattice,FullImaginary,T_Pl,T_G0,T_Gl}
    Γ = inv(G.Cw * G.insertion_trial)
    if IsLattice
        Nsub = size(G.G0)[end]
        inv(G.Cw * G.insertion_trial) * Nsub
    else
        @assert allequal(size(G.G0))
        inv(G.Cw * G.insertion_trial) * size(G.G0)[end]
    end
end

function normalize_density_matrix(G::GreenFuncBin)
    return _norm_coeff(G) .* real(G.G0)
end
# function normalize_density_matrix(Dmat::Array{Complex{Int},4}, ishardcore=true)
#     @assert size(Dmat, 3) == size(Dmat, 4)
#     Nsub = size(Dmat, 3)
#     d = [Dmat[1, 1, i, i] for i ∈ 1:Nsub]
#     D_NormCoeff::f64 = mean(imag, d) + (ishardcore ? +1 : -1) * mean(real, d)
#     Dmat_n = let D = Dmat ./ D_NormCoeff
#         for i ∈ 1:Nsub
#             D[1, 1, i, i] = real(D[1, 1, i, i])
#         end
#         Dp = real.(D) .+ imag.(D)
#         Dp
#     end
#     return Dmat_n
# end

import Base.merge
function merge(ms::Array{GreenFuncBin{IsLattice,FullImaginary,T_Pl,T_G0,T_Gl}}) where {IsLattice,FullImaginary,T_Pl,T_G0,T_Gl}
        @assert allequal(m -> m.β, ms)
        β = ms[1].β
        @assert allequal(m -> m.Cw, ms)
        Cw = ms[1].Cw
        @assert allequal(m -> size(m.G0), ms)
        _Pl = deepcopy(ms[1]._Pl)
        G0 = sum(m.G0 for m ∈ ms)
        insertion_trial = sum(m.insertion_trial for m ∈ ms)
        Gl = if FullImaginary == nothing
            @assert T_Gl == Nothing
            nothing
        else
            sum(m.Gl for m ∈ ms)
        end
        return GreenFuncBin{IsLattice,FullImaginary,T_Pl,T_G0,T_Gl}(
            β, Cw, _Pl, G0, Gl, insertion_trial
        )
end