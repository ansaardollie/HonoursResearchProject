export PopulationInfo, AgeYearRange, ModelRanges, Optional

PopulationInfo = @NamedTuple{location::String, sex::Sex}

# AgeYearRange = Union{
#     NamedTuple{(:years, :ages),Tuple{DataRange,DataRange}},
#     NamedTuple{(:years, :ages),Tuple{AbstractVector{Int},AbstractVector{Int}}},
#     NamedTuple{()}
# }

AgeYearRange = NamedTuple{(:years, :ages),Tuple{DataRange,DataRange}}

function Base.show(io::IO, ar::AgeYearRange)
    print(io, "(\n\tyears = ", repr(ar.years), ",\n\tages = ", repr(ar.ages), "\n)")
end


Optional{T} = Union{T,Nothing}