export MortalityModel, fitrange!

mutable struct MortalityModel
    population::PopulationInfo
    ranges::Stratified{AgeYearRange}
    exposures::Stratified{AgePeriodData{Float64}}
    deaths::Stratified{AgePeriodData{Float64}}
    approximatedeaths::Stratified{AgePeriodData{Float64}}
    rates::Stratified{AgePeriodData{Float64}}
    logrates::Stratified{AgePeriodData{Float64}}
    expectedlifetimes::Stratified{AgePeriodData{Float64}}
    data::Stratified{DataFrame}
    variant::Variant
    calculationmode::CalculationMode
    parameters::ModelParameters
end

function MortalityModel(;
    population::PopulationInfo,
    ranges::ModelRanges,
    data::DataFrame,
    calculation_mode::CalculationMode=CM_JULIA,
    variant::Variant=VARIANT_LC,
    kwargs...
)

    full_data = deepcopy(data)
    fit_data = subset(data, ranges.fit.years, ranges.fit.ages)
    forecast_data = subset(data, ranges.forecast.years, ranges.forecast.ages)

    strat_df = Stratified{DataFrame}(full_data, fit_data, forecast_data, "Raw Data Frames")

    df_columns = [:Exposures, :Rates, :LogRates, :Deaths, :ApproximateDeaths, :LifeExpectancies]

    ds = Dict()

    for column in df_columns
        datasets = Dict()
        for strata_sym in eachindex(strat_df)
            df_subset = strat_df[strata_sym]
            data_matrix = extract(df_subset, column)::Matrix{Float64}
            strata_range = ranges[strata_sym]

            ageyear_dataset = AgePeriodData(column, data_matrix, strata_range.ages, strata_range.years)
            datasets[strata_sym] = ageyear_dataset
        end

        stratified_data = Stratified{AgePeriodData{Float64}}(
            datasets[:full],
            datasets[:fit],
            datasets[:forecast],
            datasets[:full].label
        )

        ds[column] = stratified_data
    end

    mp = ModelParameters(ranges.fit)

    return MortalityModel(
        population,
        ranges,
        ds[:Exposures],
        ds[:Deaths],
        ds[:ApproximateDeaths],
        ds[:Rates],
        ds[:LogRates],
        ds[:LifeExpectancies],
        strat_df,
        variant,
        calculation_mode,
        mp
    )

end

InitialRangeSelection = Union{Nothing,NamedTuple{(:years, :ages),Tuple{Any,Any}}}

function MortalityModel(
    country::AbstractString,
    sex::Sex;
    remove_missing::Bool=true,
    fityears::Optional{AbstractVector{Int}}=nothing,
    fitages::Optional{AbstractVector{Int}}=nothing,
    forecastyears::Optional{AbstractVector{Int}}=nothing,
    forecastages::Optional{AbstractVector{Int}}=nothing,
    calculation_mode::CalculationMode=CM_JULIA,
    variant::Variant=VARIANT_LC,
    kwargs...
)
    cd = pwd()
    dir = "$cd/Raw Mortality Data/$country"
    exposure_df = readfile("$dir/exposures.txt"; remove_missing=remove_missing)
    deaths_df = readfile("$dir/deaths.txt"; remove_missing=remove_missing)

    popinfo::PopulationInfo = (location=country, sex=sex)

    lifetablepath = sexmatch(sex, "$dir/ltb.txt", "$dir/ltf.txt", "$dir/ltm.txt")

    lifetable_df = readfile(lifetablepath; remove_missing=remove_missing)

    years = exposure_df.Year
    ages = exposure_df.Age
    exposures = sexmatch(sex, exposure_df.Total, exposure_df.Female, exposure_df.Male)
    deaths = sexmatch(sex, deaths_df.Total, deaths_df.Female, deaths_df.Male)
    rates = lifetable_df.mx
    logrates = log.(rates)
    approxdeaths = exposures .* rates
    les = lifetable_df.ex

    df = DataFrame(
        :Year => years,
        :Age => ages,
        :Exposures => exposures,
        :Rates => rates,
        :LogRates => logrates,
        :Deaths => deaths,
        :ApproximateDeaths => approxdeaths,
        :LifeExpectancies => les
    )

    fitfor = (years=fityears, ages=fitages)
    forecastfor = (years=forecastyears, ages=forecastages)

    full_years = unique(df.Year)
    full_ages = unique(df.Age)
    fit_years = isnothing(fitfor) || isnothing(fitfor.years) ? full_years : collect(fitfor.years)
    fit_ages = isnothing(fitfor) || isnothing(fitfor.ages) ? full_ages : collect(fitfor.ages)
    fc_years = isnothing(forecastfor) || isnothing(forecastfor.years) ? full_years : collect(forecastfor.years)
    fc_ages = isnothing(forecastfor) || isnothing(forecastfor.ages) ? full_ages : collect(forecastfor.ages)

    modeldims::ModelRanges = Stratified{AgeYearRange}(
        ageyear(full_years, full_ages),
        ageyear(fit_years, fit_ages),
        ageyear(fc_years, fc_ages),
        "Data Ranges"
    )


    return MortalityModel(;
        population=popinfo,
        ranges=modeldims,
        data=df,
        calculation_mode=calculation_mode,
        variant=variant,
        kwargs...
    )

end

SecondaryRangeSelection = Union{Nothing,DataRange,AbstractVector{Int}}


function Base.show(io::IO, t::MIME"text/plain", model::MortalityModel)

    ds = displaysize(io)
    width = ds[2]

    line = repeat("=", width)
    println(io, line)
    println(io, "Mortality Model")
    println(io, "Country: ", model.population.location)
    println(io, "Sex: ", sexmatch(model.population.sex, "Males & Females", "Females", "Males"))
    println(io, line)
    println(io)
    println(io, "Data Ranges")
    println(io)

    d1ages = strip("$(model.ranges.full.ages)")
    d1years = strip("$(model.ranges.full.years)")
    d2ages = strip("$(model.ranges.fit.ages)")
    d2years = strip("$(model.ranges.fit.years)")
    d3ages = strip("$(model.ranges.forecast.ages)")
    d3years = strip("$(model.ranges.forecast.years)")

    dm = ["Full" d1ages d1years; "Fit" d2ages d2years; "Forecast" d3ages d3years]

    headers = [:Dataset, :Ages, :Years]
    model_config = table_config(io,
        headers=headers,
        alignment=:c,
        hlines=:all
    )


    println(io)
    pretty_table(io, dm; model_config...)
    println(line)

    println(io)
    println(io, "Parameters")
    println(io)

    Base.show(io, t, model.parameters)
    println(io, line)

end


