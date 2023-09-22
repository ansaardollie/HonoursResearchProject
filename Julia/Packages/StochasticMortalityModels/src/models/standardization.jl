

AgeRanges = Vector{Pair{String,Vector{Int}}}
AgeWeights = Vector{Pair{String,Vector{Float64}}}
StandardAgeRanges::AgeRanges = [
    "Infant" => [0],
    "Toddler" => [1, 2, 3, 4, 5],
    "Child" => collect(6:18),
    "Young Adult" => collect(18:25),
    "Adult" => collect(26:59),
    "Pensioner" => collect(60:80),
    "Senesence" => collect(81:110)
]

function calculate_age_weights(model::MortalityModel, ranges::AgeRanges; min_year::Optional{Int}=nothing, max_year::Optional{Int})

    all_years = years(model, DS_TRAIN)
    all_ages = ages(model, DS_TRAIN)
    low_t = isnothing(min_year) ? all_years[begin] : min_year
    high_t = isnothing(max_year) ? all_years[end] : max_year

    low_ti = indexin(low_t, all_years)[1]
    high_ti = indexin(high_t, all_years)[1]
    ext = exposures(model)[:, low_ti:high_ti]

    age_indexes = map(ranges) do p
        age_values = p[2]
        age_idx = indexin(age_values, all_ages)
        return p[1] => age_idx
    end

    range_exposures = map(age_indexes) do p
        idx = p[2]
        exposures_for_range = ext[idx, :]
        return p[1] => exposures_for_range
    end

    range_weights = map(range_exposures) do p
        range_ext = p[2]
        time_avg_ext = vec(mean(range_ext, dims=2))
        total_ext = sum(time_avg_ext)
        weights = time_avg_ext ./ total_ext
        return p[1] => weights
    end

    return range_weights
end
