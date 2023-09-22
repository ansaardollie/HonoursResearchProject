export MortalityModel,
    sex,
    location,
    calculation,
    adjustment,
    jumpoff,
    ages,
    Ages,
    years,
    Years,
    rates,
    logrates,
    deaths,
    exposures,
    lifespans,
    Rates,
    LogRates,
    Deaths,
    Exposures,
    Lifespans,
    alphas,
    Alphas,
    betas,
    Betas,
    kappas,
    Kappas

mutable struct MortalityModel
    population::PopulationInfo
    ranges::Stratified{AgeYearRange}
    dataframes::Stratified{DataFrame}
    data::Stratified{MortalityData}
    variant::ModelImplementation
    parameters::ModelParameters
    unadjusted_parameters::ModelParameters
end

function MortalityModel(;
    population::PopulationInfo,
    ranges::Stratified{AgeYearRange},
    data::DataFrame,
    variant::ModelImplementation=lc
)

    all_data = deepcopy(data)
    fit_data = subset(data, ranges.train.years, ranges.train.ages)
    test_data = subset(data, ranges.test.years, ranges.test.ages)

    strat_df = Stratified{DataFrame}(all_data, fit_data, test_data, "Raw Data Frames")

    df_columns = [:Rates, :Exposures, :Deaths, :Lifespans]

    ds = Dict{Symbol,MortalityData}()

    strata = [:complete, :train, :test]

    for s in strata
        years = ranges[s].years.values
        ages = ranges[s].ages.values
        data_dict = Dict{Symbol,Matrix{Float64}}()
        df_subset = strat_df[s]
        for c in df_columns
            vals = extract(df_subset, c)::Matrix{Float64}
            data_dict[c] = vals
        end
        rates = data_dict[:Rates]
        rates = fix_zero_rates!(rates)
        logrates = log.(rates)
        exposures = data_dict[:Exposures]
        deaths = data_dict[:Deaths]
        lifespans = data_dict[:Lifespans]
        md = MortalityData(MDS_OBSERVED, ages, years, rates, logrates, exposures, deaths, lifespans)
        ds[s] = md
    end

    strat_md = Stratified{MortalityData}(
        ds[:complete],
        ds[:train],
        ds[:test],
        "Mortality Data"
    )

    mp = ModelParameters(ranges.train)

    ump = ModelParameters(ranges.train)

    return MortalityModel(
        population,
        ranges,
        strat_df, strat_md, variant,
        mp,
        ump
    )

end


function MortalityModel(
    country::AbstractString,
    sex::Sex;
    remove_missing::Bool=true,
    train_years::Optional{AbstractVector{Int}}=nothing,
    train_ages::Optional{AbstractVector{Int}}=nothing,
    test_years::Optional{AbstractVector{Int}}=nothing,
    test_ages::Optional{AbstractVector{Int}}=nothing,
    calculation_mode::CalculationChoice=CC_JULIA,
    variant::ModelImplementation=lc_j
)
    cd = pwd()
    dir = "$cd/Raw Mortality Data/$country"
    if !isdir(dir)
        dir = "https://raw.githubusercontent.com/ansaardollie/Raw-Mortality-Data/main/$country"
    end
    exposure_df = readfile("$dir/exposures.txt"; remove_missing=remove_missing)
    deaths_df = readfile("$dir/deaths.txt"; remove_missing=remove_missing)

    popinfo::PopulationInfo = (location=country, sex=sex)

    lifetablepath = sexmatch(sex, "$dir/ltb.txt", "$dir/ltf.txt", "$dir/ltm.txt")

    lifetable_df = readfile(lifetablepath; remove_missing=remove_missing)

    years = exposure_df.Year
    ages = exposure_df.Age
    exposures = sexmatch(sex, exposure_df.Total, exposure_df.Female, exposure_df.Male)
    true_deaths = sexmatch(sex, deaths_df.Total, deaths_df.Female, deaths_df.Male)
    rates = lifetable_df.mx


    logrates = log.(rates)
    approx_deaths = exposures .* rates
    les = lifetable_df.ex

    df = DataFrame(
        :Year => years,
        :Age => ages,
        :Exposures => exposures,
        :Rates => rates,
        :LogRates => logrates,
        :Deaths => (calculation_mode == CC_JULIA ? true_deaths : approx_deaths),
        :Lifespans => les
    )

    fitfor = (years=train_years, ages=train_ages)
    testfor = (years=test_years, ages=test_ages)

    potential_years = sort(unique(df.Year))
    sy = potential_years[begin]
    ey = potential_years[end]
    all_years = Vector{Int}()

    for y in reverse(sy:ey)
        if y in potential_years
            push!(all_years, y)
        else
            break
        end
    end
    all_ages = sort(unique(df.Age))
    all_years = reverse(all_years)

    fit_years = isnothing(fitfor.years) ? all_years : collect(fitfor.years)
    fit_ages = isnothing(fitfor.ages) ? all_ages : collect(fitfor.ages)
    test_years = isnothing(testfor.years) ? all_years : collect(testfor.years)
    test_ages = !isnothing(testfor.ages) ? collect(testfor.ages) : !isnothing(fitfor.ages) ? fit_ages : all_ages



    all_ayr = ageyear(all_years, all_ages)
    fit_ayr = ageyear(fit_years, fit_ages)
    test_ayr = ageyear(test_years, test_ages)

    modeldims = Stratified{AgeYearRange}(
        all_ayr,
        fit_ayr,
        test_ayr,
        "Data Ranges"
    )

    df = DataFrames.subset(df, :Year => ByRow(y -> y in all_years))

    @reset variant.calculation = calculation_mode
    return MortalityModel(;
        population=popinfo,
        ranges=modeldims,
        data=df,
        variant=variant
    )

end



function Base.show(io::IO, t::MIME"text/plain", model::MortalityModel)

    ds = displaysize(io)
    width = ds[2]

    line = repeat("=", width)
    println(io, line)
    println(io, "Mortality Model")
    println(io, "Country: ", location(model))
    println(io, "Sex: ", sexmatch(sex(model), "Males & Females", "Females", "Males"))
    println(io, line)
    println(io)
    println(io, "Data Ranges")
    println(io)

    d1ages = strip("$(model.ranges.complete.ages)")
    d1years = strip("$(model.ranges.complete.years)")
    d2ages = strip("$(model.ranges.train.ages)")
    d2years = strip("$(model.ranges.train.years)")
    d3ages = strip("$(model.ranges.test.ages)")
    d3years = strip("$(model.ranges.test.years)")

    dm = ["All Data" d1ages d1years; "Fit Data" d2ages d2years; "Test Data" d3ages d3years]

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


sex(m::MortalityModel)::Sex = m.population.sex
location(m::MortalityModel)::String = m.population.location
calculation(m::MortalityModel)::CalculationChoice = m.variant.calculation
adjustment(m::MortalityModel)::AdjustmentChoice = m.variant.adjustment
jumpoff(m::MortalityModel)::JumpoffChoice = m.variant.jumpoff
ages(m::MortalityModel, strata::DataStrata=DS_TRAIN)::Vector{Int} = stratamatch(strata, m.data.complete.ages, m.data.train.ages, m.data.test.ages)
Ages(m::MortalityModel, strata::DataStrata=DS_TRAIN)::DataRange = stratamatch(strata, m.ranges.complete.ages, m.ranges.train.ages, m.ranges.test.ages)
years(m::MortalityModel, strata::DataStrata=DS_TRAIN)::Vector{Int} = stratamatch(strata, m.data.complete.years, m.data.train.years, m.data.test.years)
Years(m::MortalityModel, strata::DataStrata=DS_TRAIN)::DataRange = stratamatch(strata, m.ranges.complete.years, m.ranges.train.years, m.ranges.test.years)
rates(m::MortalityModel, strata::DataStrata=DS_TRAIN)::Matrix{Float64} = stratamatch(strata, m.data.complete.rates, m.data.train.rates, m.data.test.rates)

logrates(m::MortalityModel, strata::DataStrata=DS_TRAIN)::Matrix{Float64} = stratamatch(strata, m.data.complete.logrates, m.data.train.logrates, m.data.test.logrates)

exposures(m::MortalityModel, strata::DataStrata=DS_TRAIN)::Matrix{Float64} = stratamatch(strata, m.data.complete.exposures, m.data.train.exposures, m.data.test.exposures)

deaths(m::MortalityModel, strata::DataStrata=DS_TRAIN)::Matrix{Float64} = stratamatch(strata, m.data.complete.deaths, m.data.train.deaths, m.data.test.deaths)

lifespans(m::MortalityModel, strata::DataStrata=DS_TRAIN)::Matrix{Float64} = stratamatch(strata, m.data.complete.lifespans, m.data.train.lifespans, m.data.test.lifespans)

Rates(m::MortalityModel, strata::DataStrata=DS_TRAIN)::AgePeriodData{Float64} = stratamatch(strata, Rates(m.data.complete), Rates(m.data.train), Rates(m.data.test))

LogRates(m::MortalityModel, strata::DataStrata=DS_TRAIN)::AgePeriodData{Float64} = stratamatch(strata, LogRates(m.data.complete), LogRates(m.data.train), LogRates(m.data.test))

Exposures(m::MortalityModel, strata::DataStrata=DS_TRAIN)::AgePeriodData{Float64} = stratamatch(strata, Exposures(m.data.complete), Exposures(m.data.train), Exposures(m.data.test))

Deaths(m::MortalityModel, strata::DataStrata=DS_TRAIN)::AgePeriodData{Float64} = stratamatch(strata, Deaths(m.data.complete), Deaths(m.data.train), Deaths(m.data.test))

Lifespans(m::MortalityModel, strata::DataStrata=DS_TRAIN)::AgePeriodData{Float64} = stratamatch(strata, Lifespans(m.data.complete), Lifespans(m.data.train), Lifespans(m.data.test))

alphas(m::MortalityModel, version::ParameterVersion=PV_ADJUSTED)::Vector{Float64} = version == PV_ADJUSTED ? m.parameters.alphas.values : m.unadjusted_parameters.alphas.values
Alphas(m::MortalityModel, version::ParameterVersion=PV_ADJUSTED)::ParameterSet{Float64} = version == PV_ADJUSTED ? m.parameters.alphas : m.unadjusted_parameters.alphas
betas(m::MortalityModel, version::ParameterVersion=PV_ADJUSTED)::Vector{Float64} = version == PV_ADJUSTED ? m.parameters.betas.values : m.unadjusted_parameters.betas.values
Betas(m::MortalityModel, version::ParameterVersion=PV_ADJUSTED)::ParameterSet{Float64} = version == PV_ADJUSTED ? m.parameters.betas : m.unadjusted_parameters.betas
kappas(m::MortalityModel, version::ParameterVersion=PV_ADJUSTED)::Vector{Float64} = version == PV_ADJUSTED ? m.parameters.kappas.values : m.unadjusted_parameters.kappas.values
Kappas(m::MortalityModel, version::ParameterVersion=PV_ADJUSTED)::ParameterSet{Float64} = version == PV_ADJUSTED ? m.parameters.kappas : m.unadjusted_parameters.kappas

