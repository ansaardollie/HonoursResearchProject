export DataRange, ageyear, agerange, yearrange

struct DataRange <: AbstractUnitRange{Int}
    values::Vector{Int}
    start::Int
    stop::Int
    step::Union{Int,Nothing,Dict{Int,Int}}
    type::DataDimension
    string_repr::String
    indexer::Dict{Int,Int}
end


function inconsistent_range(inputs::Vector{Int})
    output = IOBuffer()
    outvals::Vector{String} = []
    years = deepcopy(inputs)
    startitem = popfirst!(years)
    while length(years) > 0
        curitem = popfirst!(years)
        enditem = curitem
        curdif = curitem - startitem
        lastdif = curitem - startitem
        while curdif == lastdif && length(years) > 0
            nextitem = popfirst!(years)
            curdif = nextitem - curitem
            if curdif != lastdif || length(years) == 0
                enditem = length(years) == 0 ? nextitem : curitem
            end
            curitem = nextitem
        end
        rangeout = lastdif == 1 ? "[$startitem to $enditem (1y intervals)]" : "[$startitem to $enditem ($(lastdif)y intervals)]"
        push!(outvals, rangeout)
        startitem = curitem
    end

    println(output)
    for r in outvals
        println(output, r)
    end

    return String(take!(output))

end


function DataRange(vals::Vector{Int}, type::DataDimension)
    start = vals[1]
    stop = vals[end]
    if length(vals) > 1
        step_sizes = vals[2:end] - vals[1:end-1]
        step = allequal(step_sizes) ? step_sizes[1] : nothing
        if allequal(step_sizes)
            step = step_sizes[1]
            repr = step > 1 ? "[$start to $stop ($(step) invervals)]" : "[$start to $stop]"
        else
            repr = inconsistent_range(vals)
        end
    elseif length(vals) == 1
        step = nothing
        repr = "[$(vals[1])]"
    else
        throw(ArgumentError("Need non-zero length for `vals`."))
    end
    idx = collect(eachindex(vals))
    idxr = Dict(zip(vals, idx))

    return DataRange(vals, start, stop, step, type, repr, idxr)
end


function DataRange(vals::OrdinalRange{Int,Int}, type::DataDimension)
    return DataRange(collect(vals), type)
end

function ageyear(years::AbstractVector{Int}, ages::AbstractVector{Int})
    ydr = DataRange(years, DD_YEAR)
    adr = DataRange(ages, DD_AGE)

    return (years=ydr, ages=adr)
end

function agerange(ages::AbstractVector{Int})
    return DataRange(ages, DD_AGE)
end

function yearrange(years::AbstractVector{Int})
    return DataRange(years, DD_YEAR)
end

function Base.print(io::IO, dr::DataRange, xs...)


    print(io, "$(dr.string_repr)")
end

function Base.show(io::IO, ::MIME"text/plain", dr::DataRange)
    Base.print(io, dr)
end


function Base.iterate(dr::DataRange)
    length(dr.values) < 1 && return nothing
    return (dr.values[1], 1)
end

function Base.iterate(dr::DataRange, state::Integer)
    length(dr.values) <= state && return nothing
    return (dr.values[state+1], state + 1)
end

function Base.IteratorSize(::DataRange)
    return Base.HasLength()
end

function Base.IteratorEltype(::DataRange)
    return Base.HasEltype()
end

function Base.eltype(::DataRange)
    return Int
end

function Base.length(dr::DataRange)
    return length(dr.values)
end

function Base.isdone(dr::DataRange)
    return length(dr.values) < 1
end


function Base.isdone(dr::DataRange, state::Integer)
    return length(dr.values) <= state
end


function Base.firstindex(dr::DataRange)
    return dr.values[begin]
end

function Base.lastindex(dr::DataRange)
    return dr.values[end]
end