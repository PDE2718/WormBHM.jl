####### Type alias ######
const f64::Type = Float64
const u64::Type = UInt64
const u128::Type = UInt128
const f32::Type = Float32
const i32::Type = Int32
const i8::Type = Int8
const StateType::Type = Int8
const IndexType::Type = Int32
const t_eps::f64 = 5e-16
metro(p)::Bool = rand() < p
randsign()::Int8 = rand(Bool) ? Int8(1) : Int8(-1)
function dumptoNamedTuple(x)
    keys = fieldnames(typeof(x))
    vals = [getfield(x, k) for k ∈ keys]
    return (; zip(keys, vals)...)
end
@generated function peel_last(x::NTuple{N,T}) where {N,T}
    return quote
        (Base.@ntuple $(N - 1) i -> x[i]), x[$(N)]
    end
end

"""
    @enum OpType

An enumeration representing different types of operations in the context of quantum mechanics or similar fields.

### Constants:
- `b_`: Represents the annihilation operation, with a value of -1.
- `I_`: Represents the identity operation, typically used as a dummy operation, with a value of 0.
- `b̂_`: Represents the creation operation, with a value of +1.
"""
@enum OpType::Int8 begin
    b_ = -1 # annihilation
    I_ = 0  # Identity (e.g. for dummy)
    b̂_ = +1 # creation
end

@enum WormLocation::Int8 begin
    _at_null # null worm tag
    _at_free # at free time on a world line
    _at_stop # head just before or after (near) tail
    _at_kink # head is near a kink
    _at_nbkink # head is near a neighboring kink. Usually can pass
    _at_dummy # head is before β or after 0. near the dummy element
    _at_green # head is of the same time as tail, but not at the same site.
end

"""
    struct Element

A structure representing an element in a system.

# Fields
- `t::f64`: The time of the element.
- `i::IndexType`: The site of the element itself.
- `j::IndexType`: The other site of a bond; 0 indicates the worm head or tail.
- `n_L::StateType`: The state to the left of the element.
- `n_R::StateType`: The state to the right of the element.
- `op::OpType`: The operator type associated with the element.

note that for comparing `Element`, the operator `>,<,==` are reloaded, and one must use `===` for the check for identity.

the operator `<<` is reloaded to change the time of a Element. e.g.:
```julia
    Element() << 1.0
```
will generate a default Element except the time is at 1.0
"""
struct Element
    t::f64          # time of the element
    i::IndexType    # the site of element it self.
    j::IndexType    # the other site of a bond, 0 for the worm head or tail.
    n_L::StateType  # state to the left  of the element.
    n_R::StateType  # state to the right of the element.
    op::OpType      # the operator type
end

Element()::Element = Element(0.,IndexType(0), IndexType(0), i8(0), i8(0), I_)
function Element(op::OpType, t::f64,
    i::IndexType, n_L::StateType, n_R::StateType,
    j::IndexType=IndexType(0))::Element
    return Element(t, i, j, n_L, n_R, op)
end
import Base: <<
<<(x::Element, t::f64) = Element(t, x.i, x.j, x.n_L, x.n_R, x.op)

const Wline::Type = Vector{Element}

dummy_element(i::IndexType, n0::StateType) = Element(1., i, IndexType(0), n0, n0, I_)

empty_wline(i0::Integer, n0::Integer) = Element[dummy_element(IndexType(i0),n0)]

function display_Wline(l::Wline)
    println("Wline : i = $(l[end].i), τ/β = $(l[end].t)")
    for e ∈ l
        println("ψ $(e.n_L)→$(e.n_R), op:$(e.op), t=$(e.t)")
    end
    println(" ")
end
import Base: isless, isequal, ==
isless(a::Element, b::Element)::Bool = isless(a.t, b.t)
isless(a::Element, t::f64)::Bool = isless(a.t, t)
isless(t::f64, b::Element)::Bool = isless(t, b.t)
==(a::Element, b::Element)::Bool = ==(a.t, b.t)
==(a::Element, t::f64)::Bool = ==(a.t, t)
==(t::f64, b::Element)::Bool = ==(t, b.t)

# function vindex(l::Wline, t::f64)::Int
#     return @inbounds searchsortedfirst(l, t)
# end
function vindex(l::Wline, t::f64)::Int
    ESIZE = sizeof(Element)
    left = 0
    right = lastindex(l) - 1
    l0 = pointer(l) |> Ptr{Float64}
    if right == -1
        return 1
    end
    if unsafe_load(muladd(ESIZE, right, l0)) < t
        return right + 2
    end
    while left < right
        mid = (left + right) >>> 1
        if unsafe_load(muladd(ESIZE, mid, l0)) < t
            left = mid + 1
        else
            right = mid
        end
    end
    return left + 1
end

function slice_at(l::Wline, t::f64)::StateType
    return l[vindex(l, t)].n_L
end

function close_to_any(l::Wline, t::f64)::Bool
    return vindex(l, t + t_eps) ≠ vindex(l, t - t_eps)
end

function check_wl(l::Wline)::Bool
    for i ∈ (1:lastindex(l)-1)
        if !(l[i].op ≠ I_ && l[i].n_R == l[i+1].n_L && l[i].t < l[i+1].t)
            return false
        end
    end
    return l[end].op == I_ && (l[end].n_L == l[end].n_R == l[1].n_L)
end

# vertex id around t. if at = true, then exlude the vertex at t (it must exist).

"""
    vindex_around(l::Wline, t::f64, at::Bool=false) :: NTuple{2, Int}

Returns a tuple containing the indices of the elements surrounding a specified time.

# Arguments
- `l::Wline`: The `Wline` object from which to retrieve the indices.
- `t::f64`: The value for which the indices are to be found.
- `at::Bool`: when `` 

# Returns
A tuple of two integers representing the indices of the elements to the left and right of the found index.

"""
@inline function vindex_around(l::Wline, t::f64, at::Bool=false)::NTuple{2,Int}
    len::Int = lastindex(l)
    i::Int = vindex(l, t)
    i_L::Int = mod1(i - 1, len)
    i_R::Int = mod1(i + at, len)
    return (i_L,i_R)
end
@inline function element_around(l::Wline, t::f64, D::Integer, at::Bool=false)::Element
    len::Int = lastindex(l)
    i::Int = vindex(l, t)
    if D > 0
        i = mod1(i + at, len)
    else
        i = mod1(i - 1, len)
    end
    return @inbounds l[i]
end

# function sizehint_wl!(l)
#     L0 = length(l)
#     if l.ref.mem.length > L0 + 8
#         sizehint!(l, L0 + 4; shrink=true)
#     # println("triggered")
#     end
#     return nothing
# end
function sizehint_wl!(l)
    tol = 2length(l)
    if l.ref.mem.length > tol > 64
        sizehint!(l, tol; shrink=true)
    end
    return nothing
end

"""
    Wsheet{N}

The imaginary time configuration of dimensions `N`, consists of many world lines.

### Fields
- `β::f64`: The inverse temperature.
- `wl::Array{Wline,N}`: An array of `Wline` (which is the alias for `Vector{Element}`), where `N` specifies the number of dimensions of the array. Here the time of each world line is scaled to `(0,1]`.
"""
struct Wsheet{N}
    β::f64
    wl::Array{Wline,N}
end
function Wsheet(β::f64, ψ0::AbstractArray{<:Integer,N})::Wsheet where {N}
    @assert β > 0
    return Wsheet(β,[empty_wline(i, ψ0[i]) for i ∈ LinearIndices(ψ0)])
end
import Base: getindex, eachindex, size, LinearIndices, CartesianIndices
function getindex(X::Wsheet{N}, i::Integer)::Wline where {N}
    return getindex(X.wl, i)
end
function getindex(X::Wsheet{N}, inds...)::Wline where {N}
    return getindex(X.wl, inds...)
end
function eachindex(X::Wsheet{N}) where {N}
    return eachindex(X.wl)
end
function LinearIndices(X::Wsheet{N}) where {N}
    return LinearIndices(X.wl)
end
function CartesianIndices(X::Wsheet{N}) where {N}
    return CartesianIndices(X.wl)
end
function size(X::Wsheet{N}) where {N}
    return size(X.wl)
end

function set_empty_wl!(X, conf)
    for (l,n0) ∈ zip(X.wl, conf)
        @assert length(l) == 1
        e = l[1]
        @assert e.op == I_ && e.j == IndexType(0)
        l[1] = Element(e.t, e.i, e.j, n0, n0, e.op)
    end
    return nothing
end
