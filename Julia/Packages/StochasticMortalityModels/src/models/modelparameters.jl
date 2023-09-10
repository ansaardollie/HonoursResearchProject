export ModelParameters, constrain!, constrain_demography!, constrain_julia!, mxt_hat, fitted_deaths

struct ModelParameters
    alphas::ParameterSet{Float64}
    betas::ParameterSet{Float64}
    kappas::ParameterSet{Float64}
    ranges::AgeYearRange
    function ModelParameters(
        alphas::ParameterSet{Float64},
        betas::ParameterSet{Float64},
        kappas::ParameterSet{Float64},
        ranges::AgeYearRange)

        ages = ranges.ages.values

        years = ranges.years.values

        sum(ages .== alphas.range.values) == length(ages) || throw(DimensionMismatch("Ages not consistent for `alphas`"))
        sum(ages .== betas.range.values) == length(ages) || throw(DimensionMismatch("Ages not consistent for `betas`"))
        sum(years .== kappas.range.values) == length(years) || throw(DimensionMismatch("Years not consistent for `kappas`"))

        alphas.range.type == DD_AGE || throw(ArgumentError("`alphas` not defined by age"))
        betas.range.type == DD_AGE || throw(ArgumentError("`betas` not defined by age"))
        kappas.range.type == DD_YEAR || throw(ArgumentError("`kappas` not defined by year"))

        return new(alphas, betas, kappas, ranges)
    end
    function ModelParameters(ranges::AgeYearRange)
        ages = ranges.ages
        years = ranges.years

        alphas = ParameterSet("α(x)", ages)
        betas = ParameterSet("β(x)", ages)
        kappas = ParameterSet("k(t)", years)


        return new(alphas, betas, kappas, ranges)
    end
end

function ModelParameters(alphas::ParameterSet{Float64}, betas::ParameterSet{Float64}, kappas::ParameterSet{Float64})
    return ModelParameters(alphas, betas, kappas, (years=kappas.range, ages=alphas.range))
end

function ModelParameters(alpha_values::Vector{Float64}, beta_values::Vector{Float64}, kappa_values::Vector{Float64}, ages::AbstractVector{Int}, years::AbstractVector{Int})

    ar = agerange(ages)
    yr = yearrange(years)
    alphas = ParameterSet("Alpha", alpha_values, ar)
    betas = ParameterSet("Beta", beta_values, ar)
    kappas = ParameterSet("Kappa", kappa_values, yr)

    return ModelParameters(alphas, betas, kappas, (years=yr, ages=ar))
end


function ModelParameters(; alphas::Vector{Float64}, betas::Vector{Float64}, kappas::Vector{Float64}, ages::Union{DataRange,Vector{Int}}, years::Union{DataRange,Vector{Int}})
    ar = ages isa DataRange ? ages : agerange(ages)
    yr = years isa DataRange ? years : yearrange(years)

    return ModelParameters(alphas, betas, kappas, ar.values, yr.values)
end



@inline function Base.getindex(mp::ModelParameters, i::Int)
    @boundscheck 1 <= i <= 3 || throw_boundserror(mp, i)
    outcome = @match i begin
        1 => mp.alphas
        2 => mp.betas
        3 => mp.kappas
    end
    return outcome
end

@inline function Base.getindex(mp::ModelParameters, i::Symbol)
    @boundscheck in(i, [:Alphas, :Betas, :Kappas, :alphas, :betas, :kappas, :alpha, :beta, :kappa, :Alpha, :Beta, :Kappa, :α, :β, :κ]) || throw_boundserror(mp, i)

    indx = @match i begin
        :Alphas || :Alpha || :alphas || :alpha || :α => 1
        :Betas || :Beta || :betas || :beta || :β => 2
        :Kappas || :Kappa || :kappas || :kappa || :κ => 3
    end

    return @inbounds mp[indx]
end

@inline function Base.setindex!(mp::ModelParameters, values::Vector{Float64}, i::Int)
    @boundscheck 1 <= i <= 3 || throw_boundserror(mp, i)

    req_length = i == 3 ? length(mp.ranges.years) : length(mp.ranges.ages)

    length(values) == req_length || throw(ArgumentError("Input vector size ($(length(values))) does not match parameter set size ($req_length)"))

    if i == 1
        @inbounds mp.alphas.values[:] = values
    elseif i == 2
        @inbounds mp.betas.values[:] = values
    else
        @inbounds mp.kappas.values[:] = values
    end

    return mp
end

@inline function Base.setindex!(mp::ModelParameters, values::Vector{Float64}, i::Symbol)
    @boundscheck in(i, [:Alphas, :Betas, :Kappas, :alphas, :betas, :kappas, :alpha, :beta, :kappa, :Alpha, :Beta, :Kappa, :α, :β, :κ]) || throw_boundserror(mp, i)

    indx = @match i begin
        :Alphas || :Alpha || :alphas || :alpha || :α => 1
        :Betas || :Beta || :betas || :beta || :β => 2
        :Kappas || :Kappa || :kappas || :kappa || :κ => 3
    end

    return @inbounds Base.setindex!(mp, values, indx)
end

@inline function Base.setindex!(mp::ModelParameters, params::ParameterSet{Float64}, i::Int)
    return Base.setindex!(mp, params.values, i)
end

@inline function Base.setindex!(mp::ModelParameters, params::ParameterSet{Float64}, i::Symbol)
    return Base.setindex!(mp, params.values, i)
end


function Base.show(io::IO, t::MIME"text/plain", mp::ModelParameters)

    ages = mp.alphas.range.values
    years = mp.kappas.range.values

    title_out = "Age Dependent Parameters"
    title_year = "Time Dependent Parameters"


    age_conf = table_config(
        io,
        headers=map(Symbol, [mp.alphas.name, mp.betas.name]),
        rows=ages,
        row_label_title="Ages",
        title=title_out,
        formatters=seperator
    )

    pretty_table(io, hcat(mp.alphas.values, mp.betas.values); age_conf...)

    year_conf = table_config(
        io,
        headers=map(Symbol, [mp.kappas.name]),
        rows=years,
        row_label_title="Years",
        formatters=seperator,
        title=title_year
    )

    println(io)
    println(io)
    pretty_table(io, hcat(mp.kappas.values); year_conf...)

end


function constrain!(mp::ModelParameters; mode::CalculationMode=CM_JULIA)
    if mode == CM_JULIA
        return constrain_julia!(mp)
    elseif mode == CM_DEMOGRAPHY
        return constrain_demography!(mp)
    else
        throw(ArgumentError("CalulcationMode `$mode` not understood"))
    end
end

function constrain_demography!(mp::ModelParameters)
    c1 = sum(mp.betas)

    mp[:beta] = mp[:beta] / c1
    mp[:kappa] = mp[:kappa] * c1

    return mp
end

function constrain_julia!(mp::ModelParameters)
    c1 = mean(mp.kappas)
    c2 = sum(mp.betas)

    mp[:alpha] = mp.alphas + c1 * mp.betas
    mp[:beta] = mp.betas / c2
    mp[:kappa] = c2 * (mp.kappas - c1)

    return mp
end

function mxt_hat(mp::ModelParameters; ages::Union{AbstractVector{Int},Colon}=:, years::Union{AbstractVector{Int},Colon}=:, log_scale::Bool=false)
    alpha = mp.alphas[ages]
    beta = mp.betas[ages]
    kappa = mp.kappas[years]

    lmxt = alpha + beta * kappa

    if log_scale
        @reset lmxt.category = MDC_LOGRATES
        @reset lmxt.source = MDS_FITTED
        @reset lmxt.label = mdc_label(lmxt.category, lmxt.source)
        return lmxt
    else
        mxt = map(l -> exp(l), lmxt)

        @reset mxt.category = MDC_RATES
        @reset mxt.source = MDS_FITTED
        @reset mxt.label = mdc_label(mxt.category, mxt.source)
        return mxt

    end
end

function fitted_deaths(mp::ModelParameters, exposures::AgePeriodData{Float64}; ages::Union{AbstractVector{Int},Colon}=:, years::Union{AbstractVector{Int},Colon}=:)
    mxt = mxt_hat(mp; ages=ages, years=years)
    Ext = exposures[ages, years]

    dxt::AgePeriodData{Float64} = mxt * Ext

    dxt = Accessors.setproperties(dxt, (category=MDC_DEATHS, source=MDS_FITTED, label=mdc_shortlabel(MDC_DEATHS, MDS_FITTED)))

    return dxt
end
