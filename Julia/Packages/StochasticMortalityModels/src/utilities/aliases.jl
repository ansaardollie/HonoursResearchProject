export PopulationInfo, AgeYearRange, ModelRanges, Optional

PopulationInfo = @NamedTuple{location::String, sex::Sex}

AgeYearRange = Union{
    NamedTuple{(:years, :ages),Tuple{DataRange,DataRange}},
    NamedTuple{(:years, :ages),Tuple{AbstractVector{Int},AbstractVector{Int}}},
    NamedTuple{()}
}

ModelRanges = Stratified{AgeYearRange}

Optional{T} = Union{T,Nothing}