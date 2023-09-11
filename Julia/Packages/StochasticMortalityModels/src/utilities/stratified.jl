export Stratified

struct Stratified{T} <: AbstractVector{T}
    all::T
    fit::T
    test::T
    label::Union{String,Nothing}
end

function Stratified(all::T, fit::T, test::T, label=nothing) where {T}
    return Stratified(all, fit, test, label)
end


function Base.getindex(s::Stratified{T}, i::Integer) where {T}
    outcome = @match i begin
        1 => s.all
        2 => s.fit
        3 => s.test
        _ => throw(BoundsError(s, i))
    end
    return outcome
end

function Base.getindex(s::Stratified{T}, i::Symbol) where {T}
    outcome = @match i begin
        :all => s[1]
        :fit => s[2]
        :test => s[3]
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
    if isdefined(s, :all)
        return (s[:all], :all)
    else
        return nothing
    end
end

function Base.iterate(s::Stratified{T}, state) where {T}
    if isdefined(s, state)
        return (s[state], state == :all ? :fit : :test)
    else
        return nothing
    end
end

function Base.eachindex(s::Stratified{T}) where {T}
    return [:all, :fit, :test]
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
    label = isnothing(ayd.label) ? "$(typeof(ayd.all))" : ayd.label

    print(IOContext(io))
    println(io, line)
    println(io)
    println(io, align_string("$(label) Stratified By Dataset", width, :c))
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("All Data", width, :l))
    println(io)
    Base.show(io, t, ayd.all)
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("Fitting Data", width, :l))
    println(io)
    Base.show(io, t, ayd.fit)
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("Testing Data", width, :l))
    println(io)
    Base.show(io, t, ayd.test)
    println(io)
    println(io, line)
end