mutable struct SimpleMeasure{Np}
    props::NTuple{Np,f64}
    const names::NTuple{Np,Symbol}
    n_measure::Int
end
import Base: empty!
function empty!(m::SimpleMeasure)
    m.props = 0.0 .* m.props
    m.n_measure = 0
    return m
end
import Base: show
function Base.show(io::IO, m::SimpleMeasure)
    println(io, "SimpleMeasure: ", NamedTuple(zip(m.names, m.props ./ m.n_measure)))
    nothing
end
import Base: getindex
function getindex(m::SimpleMeasure{Np}, o::Symbol) where {Np}
    i = findfirst(==(o), m.names)
    return isnothing(i) ? NaN : m.props[i] / m.n_measure
end
import Base: merge
function merge(ms::Array{SimpleMeasure{Np}}) where {Np}
    @assert allequal(m -> m.names, ms)
    names = ms[1].names
    m_merge = SimpleMeasure{Np}(0.0 .* ms[1].props, names, 0)
    for m ∈ ms
        m_merge.props = (m_merge.props .+ m.props)
        m_merge.n_measure += m.n_measure
    end
    return m_merge
end
