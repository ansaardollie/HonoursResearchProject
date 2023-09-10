using Base: throw_boundserror, setindex

export ParameterSet

struct ParameterSet{T<:Real} <: AbstractVector{T}
    name::String
    values::Vector{T}
    range::DataRange
    ParameterSet{T}(name::String, values::Vector{T}, range::DataRange) where {T} = new{T}(name, values, range)
end

function ParameterSet(name::String, values::Vector{T}, range::DataRange) where {T}
    return ParameterSet{T}(name, values, range)
end
function ParameterSet{S}(name::String, range::DataRange) where {S<:Real}
    n = length(range)
    vals = Vector{S}(undef, n)
    return ParameterSet(name, vals, range)
end
function ParameterSet(name::String, range::DataRange)
    n = length(range)
    vals = Vector{Float64}(undef, n)
    return ParameterSet(name, vals, range)
end

Base.length(ps::ParameterSet) = length(ps.values)
Base.size(ps::ParameterSet) = (length(ps),)

@inline function Base.getindex(ps::ParameterSet, i::Int)
    @boundscheck issubset(i, ps.range.values) || throw_boundserror(ps, i)

    val = ps.values[ps.range.indexer[i]]
    dr = ps.range.type == DD_AGE ? agerange(i:i) : yearrange(i:i)

    return ParameterSet(ps.name, [val], dr)
end

@inline function Base.getindex(ps::ParameterSet, I...)
    if I[1] isa Colon
        return ps
    end

    @boundscheck begin
        length(I) == 1 || throw(BoundsError(ps, I))
        issubset(I[1], ps.range.values) || throw(BoundsError(ps, I[1]))
    end

    idx = map(i -> ps.range.indexer[i], I[1])
    dr = ps.range.type == DD_AGE ? agerange(ps.range[idx]) : yearrange(ps.range[idx])

    return ParameterSet(ps.name, ps.values[idx], dr)
end

@inline function Base.setindex!(ps::ParameterSet{T}, v::T, i::Int) where {T<:Real}
    @boundscheck haskey(ps.range.indexer, i) || throw_boundserror(ps, i)

    ps.values[ps.range.indexer[i]] = v
end

function Base.setindex!(ps::ParameterSet{T}, V::Vector{T}, I...) where {T<:Real}
    @boundscheck begin
        length(I) == 1 || throw(BoundsError(ps, I))
        issubset(I[1], ps.range.values) || throw(BoundsError(ps, I[1]))
    end

    idx = map(i -> ps.range.indexer[i], I[1])
    ps.values[idx] = V
end


function Base.eachindex(ps::ParameterSet)
    return ps.range
end

function Base.axes(ps::ParameterSet)
    return (ps.range,)
end

function seperator(v, i, j)

    outcome = @sprintf "%.6f" v
    parts = split(outcome, letter -> letter == '.')

    lhs = reverse(parts[1])
    rhs = parts[2]

    lseperated = reverse(join(
        [(i % 3 == 0) ? "$(lhs[i]) " : "$(lhs[i])" for i in eachindex(lhs)]
    ))
    rseperated = join([(i % 3 == 0) ? "$(rhs[i]) " : "$(rhs[i])" for i in eachindex(rhs)])

    return "$lseperated.$rseperated"
end

function Base.show(io::IO, t::MIME"text/plain", ps::ParameterSet{T}) where {T}

    indices = ps.range.values
    indextype = ps.range.type == DD_AGE ? "Age" : "Year"


    title_out = "$(ps.name) Parameter Values\n\n"
    param_config = table_config(io,
        title=title_out,
        headers=map(Symbol, [ps.name]),
        rows=indices,
        row_label_title=indextype,
        formatters=seperator
    )

    pretty_table(io, ps.values; param_config...)
end


function Base.firstindex(ps::ParameterSet)
    return ps.range.start
end

function Base.lastindex(ps::ParameterSet)
    return ps.range.stop
end

function Base.iterate(ps::ParameterSet)
    # out = iterate(ps.range)
    # isnothing(out) && return nothing
    # (idx, st) = out
    length(ps.values) > 0 || return nothing
    return (ps.values[1], 1)
end

function Base.iterate(ps::ParameterSet, idx)
    length(ps.values) > idx || return nothing
    return (ps.values[idx+1], idx + 1)
end

MyDims = Union{Dims,Tuple{Vararg{DataRange,N}}} where {N}

function Base.print(io::IO, ps::ParameterSet, xs...)
    Base.print(io, "$(ps.name): $(ps.values) for $(ps.range.type == DD_AGE ? "Ages" : "Years") $(ps.range)")
    # Base.print(io, ps.name,": ")
    # Base.print(io,)
end

function Base.similar(ps::ParameterSet{T}, el::Type{S}, dim::MyDims) where {T<:Real,S<:Real}


    println("Similar called on ParameterSet")
    println("First arg type == $(typeof(ps))")
    println(ps)
    println("Input for `el` = $(el) with type `$(typeof(el))`.")
    println("Input for `dims` = $(dim) with type `$(typeof(dim))`")

    return ParameterSet{S}(ps.name, dim[1])

end

function Base.LinearIndices(ps::ParameterSet{T}) where {T<:Real}
    return ps.range.values
end

function Base.:*(ps::ParameterSet{T}, n::Number) where {T<:Real}
    vals = ps.values .* n
    return ParameterSet(ps.name, vals, ps.range)
end

function Base.:*(n::Number, ps::ParameterSet{T}) where {T<:Real}
    return Base.:*(ps, n)
end


function Base.:+(ps::ParameterSet{T}, n::Number) where {T<:Real}
    vals = ps.values .+ n
    return ParameterSet(ps.name, vals, ps.range)
end

function Base.:+(n::Number, ps::ParameterSet{T}) where {T<:Real}
    return Base.:+(ps, n)
end


function Base.:/(ps::ParameterSet{T}, n::Number) where {T<:Real}
    return ps * (1 / n)
end

function Base.:/(n::Number, ps::ParameterSet{T}) where {T<:Real}
    vals = n ./ ps.values
    return ParameterSet(ps.name, vals, ps.range)
end


function Base.:-(ps::ParameterSet{T}, n::Number) where {T<:Real}
    return ps + (-n)
end

function Base.:-(n::Number, ps::ParameterSet{T}) where {T<:Real}
    return n + (-1 * ps)
end


function Base.:+(ps1::ParameterSet{T}, ps2::ParameterSet{T}) where {T<:Real}
    length(ps1) == length(ps2) || throw(DimensionMismatch("Can't mix $(length(ps1)) with $(length(ps2))"))
    sum(ps1.range.values .== ps2.range.values) == length(ps1) || throw(DimensionMismatch("Can't mix $(length(ps1)) with $(length(ps2))"))
    n = "$(ps1.name) + $(ps2.name)"
    return ParameterSet(n, ps1.values .+ ps2.values, ps1.range)
end

function Base.:-(ps1::ParameterSet{T}, ps2::ParameterSet{T}) where {T<:Real}
    return ps1 + (-1 * ps2)
end

function Base.:*(ps1::ParameterSet{T}, ps2::ParameterSet{T}) where {T<:Real}
    if ps1.range.type == ps2.range.type
        length(ps1) == length(ps2) || throw(DimensionMismatch("Can't mix $(length(ps1)) with $(length(ps2))"))
        sum(ps1.range.values .== ps2.range.values) == length(ps1) || throw(DimensionMismatch("Can't mix $(length(ps1)) with $(length(ps2))"))

        newvals = ps1.values .* ps2.values

        n = "$(ps1.name) × $(ps2.name)"
        return ParameterSet(n, newvals, ps1.range)

    else
        yearparams = ps1.range.type == DD_YEAR ? ps1 : ps2
        ageparams = ps1.range.type == DD_AGE ? ps1 : ps2

        row = reshape(yearparams.values, 1, length(yearparams))
        col = reshape(ageparams.values, length(ageparams), 1)

        dm = col * row

        name = "$(ageparams.name)(x) × $(yearparams.name)(t)"
        return AgePeriodData(MDC_PARAMETERS, dm, ageparams.range, yearparams.range; rounding=6, label=name, source=MDS_CALCULATED)

    end
end

function Base.:+(ps::ParameterSet, apd::AgePeriodData)
    length(ps) == length(apd.periods) ||
        length(ps) == length(apd.ages) ||
        throw(DimensionMismatch("Size of AgePeriodData $(size(apd)) does not match input $(length(ps))"))

    direction = ps.range.type

    range_to_check = direction == DD_AGE ? apd.ages : apd.periods

    sum(ps.range.values .== range_to_check.values) == length(ps) || throw(DimensionMismatch("$(direction == DD_AGE ? "Ages" : "Years") dont match"))

    rs_size = direction == DD_YEAR ? (1, length(ps)) : (length(ps), 1)
    dm = apd.data
    vals = reshape(ps.values, rs_size...)

    result = vals .+ dm
    n = "$(ps.name)($(direction == DD_AGE ? "x" : "t")) + $(apd.label)"

    result = AgePeriodData(MDC_PARAMETERS, result, apd.ages, apd.periods; rounding=6, label=n, source=MDS_CALCULATED)
end

function Base.:+(ayd::AgePeriodData, ps::ParameterSet)
    return ps + ayd
end

