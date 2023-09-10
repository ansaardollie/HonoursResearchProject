export Stratified

struct Stratified{T} <: AbstractVector{T}
    full::T
    fit::T
    forecast::T
    label::Union{String,Nothing}
end

function Stratified(full::T, fit::T, forcast::T, label=nothing) where {T}
    return Stratified(full, fit, forecast, label)
end

function Base.getindex(s::Stratified{T}, i::Integer) where {T}
    outcome = @match i begin
        1 => s.full
        2 => s.fit
        3 => s.forecast
        _ => throw(BoundsError(s, i))
    end
    return outcome
end

function Base.getindex(s::Stratified{T}, i::Symbol) where {T}
    outcome = @match i begin
        :full => s[1]
        :fit => s[2]
        :forecast || :fc => s[3]
        _ => throw(BoundsError(s, i))
    end
    return outcome
end

function Base.firstindex(s::Stratified{T}) where {T}
    return 1
end


function Base.lastindex(s::Stratified{T}) where {T}
    return 3
end


function Base.iterate(s::Stratified{T}) where {T}
    if isdefined(s, :full)
        return (s[:full], :full)
    else
        return nothing
    end
end

function Base.iterate(s::Stratified{T}, state) where {T}
    if isdefined(s, state)
        return (s[state], state == :full ? :fit : :forecast)
    else
        return nothing
    end
end

function Base.eachindex(s::Stratified{T}) where {T}
    return [:full, :fit, :forecast]
end


function Base.size(s::Stratified{T}) where {T}
    return (3)
end


function IndexStyle(s::Stratified{T}) where {T}
    return IndexCartesian()
end


function Base.show(io::IO, t::MIME"text/plain", ayd::Stratified{T}) where {T}

    ds = displaysize(io)
    width = ds[2]

    line = repeat("=", width)
    label = isnothing(ayd.label) ? "$(typeof(ayd.full))" : ayd.label

    print(IOContext(io))
    println(io, line)
    println(io)
    println(io, align_string("Stratified data of $(label)", width, :c))
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("Full period.", width, :l))
    println(io)
    Base.show(io, t, ayd.full)
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("Fitting period", width, :l))
    println(io)
    Base.show(io, t, ayd.fit)
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("Forcasting period", width, :l))
    println(io)
    Base.show(io, t, ayd.forecast)
    println(io)
    println(io, line)
end