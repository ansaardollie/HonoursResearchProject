export Stratified

struct Stratified{T} <: AbstractVector{T}
    complete::T
    train::T
    test::T
    label::Union{String,Nothing}
end

function Stratified(all::T, train::T, test::T, label=nothing) where {T}
    return Stratified(all, train, test, label)
end


function Base.getindex(s::Stratified{T}, i::Integer) where {T}
    outcome = @match i begin
        1 => s.complete
        2 => s.train
        3 => s.test
        _ => throw(BoundsError(s, i))
    end
    return outcome
end

function Base.getindex(s::Stratified{T}, i::Symbol) where {T}
    outcome = @match i begin
        :complete => s[1]
        :train => s[2]
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
    if isdefined(s, :complete)
        return (s[:complete], :complete)
    else
        return nothing
    end
end

function Base.iterate(s::Stratified{T}, state) where {T}
    if isdefined(s, state)
        return (s[state], state == :all ? :train : :test)
    else
        return nothing
    end
end

function Base.eachindex(s::Stratified{T}) where {T}
    return [:complete, :train, :test]
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
    label = isnothing(ayd.label) ? "$(typeof(ayd.complete))" : ayd.label

    print(IOContext(io))
    println(io, line)
    println(io)
    println(io, align_string("$(label) Stratified Datasets", width, :c))
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("Complete Dataset", width, :l))
    println(io)
    Base.show(io, t, ayd.complete)
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("Training Dataset", width, :l))
    println(io)
    Base.show(io, t, ayd.train)
    println(io)
    println(io, line)
    println(io)
    println(io, align_string("Testing Dataset", width, :l))
    println(io)
    Base.show(io, t, ayd.test)
    println(io)
    println(io, line)
end