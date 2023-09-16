
struct MortalityData
    source::MortalityDataSource
    ages::Vector{Int}
    years::Vector{Int}
    rates::Matrix{Float64}
    logrates::Matrix{Float64}
    exposures::Matrix{Float64}
    deaths::Matrix{Float64}
    lifespans::Matrix{Float64}
end


ages(md)::Vector{Int} = md.ages
years(md)::Vector{Int} = md.years
Ages(md)::DataRange = agerange(md.ages)
Years(md)::DataRange = yearrange(md.years)

Rates(md)::AgePeriodData{Float64} = AgePeriodData(MDC_RATES, md.rates, Ages(md), Years(md); source=md.source)
LogRates(md)::AgePeriodData{Float64} = AgePeriodData(MDC_LOGRATES, md.logrates, Ages(md), Years(md); source=md.source)
Exposures(md)::AgePeriodData{Float64} = AgePeriodData(MDC_EXPOSURES, md.exposures, Ages(md), Years(md); source=md.source)
Deaths(md)::AgePeriodData{Float64} = AgePeriodData(MDC_DEATHS, md.deaths, Ages(md), Years(md); source=md.source)
Lifespans(md)::AgePeriodData{Float64} = AgePeriodData(MDC_LIFETIMES, md.lifetimes, Ages(md), Years(md); source=md.source)
