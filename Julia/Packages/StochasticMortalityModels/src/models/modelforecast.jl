export ForecastedData, ModelForecasts, forecast

struct ForecastedData{T}
    forecast::T
    lower::Optional{T}
    upper::Optional{T}
    label::AbstractString
    cl::Optional{Float64}
    uncertainty::Optional{UncertaintyMode}
end



struct ModelForecasts
    kappas::ForecastedData{ParameterSet{Float64}}
    rates::ForecastedData{AgePeriodData{Float64}}
    lifespans::ForecastedData{AgePeriodData{Float64}}
    deaths::ForecastedData{AgePeriodData{Float64}}
    uncertainty::Optional{UncertaintyMode}
    function ModelForecasts(
        kappas::Vector{ParameterSet{Float64}},
        rates::Vector{AgePeriodData{Float64}},
        ls::Vector{AgePeriodData{Float64}},
        deaths::Vector{AgePeriodData{Float64}},
        cl::Optional{Float64}=nothing;
        uncertainty::Optional{UncertaintyMode}=nothing
    )
        k = ForecastedData(kappas[1], kappas[2], kappas[3], "κ(t)", cl, uncertainty)
        mxt = ForecastedData(rates[1], rates[2], rates[3], "m(x,t)", cl, uncertainty)
        elt = ForecastedData(ls[1], ls[2], ls[3], "e(x,t)", cl, uncertainty)
        dxt = ForecastedData(deaths[1], deaths[2], deaths[3], "d(x,t)", cl, uncertainty)

        return new(k, mxt, elt, dxt, uncertainty)
    end
end


function forecast(
    model::MortalityModel;
    confidence_level::Float64=0.95,
    uncertainty::UncertaintyMode=UM_INNOVATION_ONLY
)
    test_years = years(model, DS_TRAIN)

    κ = kappas(model)

    use_actual = jumpoff(model) == JR_ACTUAL

    jumpoff_rates = use_actual ?
                    rates(model)[:, end] :
                    vec(mxt_hat(model.parameters, year_subset=[test_years[end]]).values)


    log_jr = log.(jumpoff_rates)

    kt = κ .- κ[end]

    drift_est = (kt[end] - kt[begin]) / (test_years[end] - test_years[begin])

    difs = kt[2:end] - kt[1:end-1]

    rw_var_est = sum((difs .- drift_est) .^ 2) / (length(difs) - 1)

    horizon_agerange = Ages(model, DS_TEST) #model.ranges.test.ages
    ages = horizon_agerange.values
    horizon_yearrange = Years(model, DS_TEST)#model.ranges.test.years
    horizon = horizon_yearrange.values

    se_σ = sqrt(rw_var_est)
    se_μ = sqrt(rw_var_est / length(difs))

    h = length(horizon)
    x = collect(1:h)

    kt_se = uncertainty == UM_INNOVATION_ONLY ? sqrt.(x .* se_σ^2) : sqrt.(x .* se_σ^2 + (x .* se_μ) .^ 2)

    α = 1 - confidence_level
    z = quantile(Normal(), 1 - (α / 2))

    kt_forecast = kt[end] .+ (x .* drift_est)
    kt_flb = kt_forecast .- (z .* kt_se)
    kt_fub = kt_forecast .+ (z .* kt_se)

    alphas = ParameterSet("Jump Off α(x)", log_jr, model.parameters.alphas.range)
    betas = model.parameters.betas



    kappa_fc = ParameterSet("κ(t) [Forecast]", kt_forecast, horizon_yearrange)
    kappa_flb = ParameterSet("κ(t) [$(round(confidence_level*100,digits=0))% LB]", kt_flb, horizon_yearrange)
    kappa_fub = ParameterSet("κ(t) [$(round(confidence_level*100,digits=0))% UB]", kt_fub, horizon_yearrange)

    forecasted_params = ModelParameters(deepcopy(alphas), deepcopy(betas), deepcopy(kappa_fc))

    flb_params = ModelParameters(deepcopy(alphas), deepcopy(betas), deepcopy(kappa_flb))

    fub_params = ModelParameters(deepcopy(alphas), deepcopy(betas), deepcopy(kappa_fub))


    mxt_fc = mxt_hat(forecasted_params)
    mxt_flb = mxt_hat(flb_params)
    mxt_fub = mxt_hat(fub_params)

    for i in eachindex(mxt_flb.values)
        lb = mxt_flb.values[i]
        ub = mxt_fub.values[i]

        if ub < lb
            mxt_flb.values[i] = ub
            mxt_fub.values[i] = lb
        end
    end



    @reset mxt_fc.source = MDS_PREDICTED
    @reset mxt_flb.source = MDS_PREDICTED
    @reset mxt_fub.source = MDS_PREDICTED
    rl = mdc_shortlabel(MDC_RATES, MDS_OBSERVED)
    @reset mxt_fc.label = "Forecasted $rl"
    @reset mxt_flb.label = "Forecasted $(round(confidence_level*100,digits=0))% LB $rl"
    @reset mxt_fub.label = "Forecasted $(round(confidence_level*100,digits=0))% LB $rl"

    cm = calculation(model)
    le_fc_dm = lexpectancies(mxt_fc.values, ages, horizon, gender=model.population.sex, at_age=ages, mode=cm)
    le_flb_dm = lexpectancies(mxt_fub.values, ages, horizon, gender=model.population.sex, at_age=ages, mode=cm)
    le_fub_dm = lexpectancies(mxt_flb.values, ages, horizon, gender=model.population.sex, at_age=ages, mode=cm)

    ext = exposures(model, DS_TEST)
    dxt_fc_dm = ext .* mxt_fc.values
    dxt_flb_dm = ext .* mxt_flb.values
    dxt_fub_dm = ext .* mxt_fub.values

    rl = mdc_shortlabel(MDC_DEATHS, MDS_OBSERVED)
    dxt_fc = AgePeriodData(MDC_DEATHS, MDS_PREDICTED, "Forecasted $rl", dxt_fc_dm, horizon_agerange, horizon_yearrange, 3, false, nothing)

    dxt_flb = AgePeriodData(MDC_DEATHS, MDS_PREDICTED, "Forecasted $(round(confidence_level*100,digits=0))% LB $rl", dxt_flb_dm, horizon_agerange, horizon_yearrange, 3, false, nothing)

    dxt_fub = AgePeriodData(MDC_DEATHS, MDS_PREDICTED, "Forecasted $(round(confidence_level*100,digits=0))% UB $rl", dxt_fub_dm, horizon_agerange, horizon_yearrange, 3, false, nothing)


    rl = mdc_shortlabel(MDC_LIFESPANS, MDS_OBSERVED)
    le_fc = AgePeriodData(MDC_LIFESPANS, MDS_PREDICTED, "Forecasted $rl", le_fc_dm, horizon_agerange, horizon_yearrange, 3, false, nothing)

    le_flb = AgePeriodData(MDC_LIFESPANS, MDS_PREDICTED, "Forecasted $(round(confidence_level*100,digits=0))% LB $rl", le_flb_dm, horizon_agerange, horizon_yearrange, 3, false, nothing)

    le_fub = AgePeriodData(MDC_LIFESPANS, MDS_PREDICTED, "Forecasted $(round(confidence_level*100,digits=0))% UB $rl", le_fub_dm, horizon_agerange, horizon_yearrange, 3, false, nothing)

    kappa = [kappa_fc, kappa_flb, kappa_fub]
    mxt = [mxt_fc, mxt_flb, mxt_fub]
    elt = [le_fc, le_flb, le_fub]
    dxt = [dxt_fc, dxt_flb, dxt_fub]

    return ModelForecasts(kappa, mxt, elt, dxt, confidence_level; uncertainty=uncertainty)

end

function Base.show(io::IO, t::MIME"text/plain", fd::ForecastedData{ParameterSet{Float64}})

    if !isnothing(fd.lower) && isnothing(fd.upper)
        dm = hcat(fd.forecast.values)
        headers = [:Forecast]
    elseif isnothing(fd.lower) && !isnothing(fd.upper)
        dm = hcat(fd.forecast.values, fd.upper.values)
        headers = [:Forcast, Symbol("Upper Bound")]
    elseif isnothing(fd.upper) && !isnothing(fd.lower)
        dm = hcat(fd.lower.values, fd.forecast.values)
        headers = [Symbol("Lower Bound"), :Forcast]
    else
        dm = hcat(fd.lower.values, fd.forecast.values, fd.upper.values)
        headers = [Symbol("Lower Bound"), Symbol("Forecast"), Symbol("Upper Bound")]
    end

    rows = fd.forecast.range.values
    rt = fd.forecast.range.type == DD_AGE ? "Ages" : "Years"

    title_io = IOBuffer()
    println(title_io, "Forecasted Data: ", fd.label)
    if !isnothing(fd.cl)
        println(title_io, "Prediction Interval Level: ", round(fd.cl * 100, digits=0), "%")
    end
    if !isnothing(fd.uncertainty)
        println(title_io, "Uncertainty Accounted For: ", fd.uncertainty == UM_INNOVATION_ONLY ? "Random Walk Innovation Only" : "Random Walk Innovation + Drift")
    end

    title = String(take!(title_io))
    fd_config = table_config(io, title=title, rows=rows, row_label_title=rt, headers=headers, formatters=gen_seperator(3))

    pretty_table(io, dm; fd_config...)


end



function Base.show(io::IO, t::MIME"text/plain", fd::ForecastedData{AgePeriodData{Float64}})

    fc_vals = filter(x -> x != 0.0, fd.forecast.values)
    lb_vals = filter(x -> x != 0.0, fd.lower.values)
    ub_vals = filter(x -> x != 0.0, fd.upper.values)
    fc_exp = Int(floor(minimum(log10.(abs.(fc_vals)))))
    ub_exp = isnothing(fd.upper) ? 0 : Int(floor(minimum(log10.(abs.(ub_vals)))))
    lb_exp = isnothing(fd.upper) ? 0 : Int(floor(minimum(log10.(abs.(lb_vals)))))

    rescale = fc_exp < -3 || ub_exp < -3 || lb_exp < -3
    scale_exp = minimum([fc_exp, ub_exp, lb_exp])
    fc_data = rescale ? Int.(round.(fd.forecast.values ./ 10.0^(scale_exp), digits=0)) : fd.forecast.values
    ub_data = isnothing(fd.upper) ? nothing : rescale ? Int.(round.(fd.upper.values ./ 10.0^(scale_exp), digits=0)) : fd.upper.values
    lb_data = isnothing(fd.lower) ? nothing : rescale ? Int.(round.(fd.lower.values ./ 10.0^(scale_exp), digits=0)) : fd.lower.values

    dm = Matrix{String}(undef, size(fc_data)...)
    sep = gen_seperator(rescale ? 0 : 2)
    for i in eachindex(fc_data)
        pe = sep(fc_data[i], 0, 0)
        if isnothing(ub_data) && isnothing(lb_data)
            dm[i] = sep(pe, 0, 0)
        elseif !isnothing(ub_data) && isnothing(lb_data)
            ub = sep(ub_data[i], 0, 0)
            cell = IOBuffer()
            print(cell, pe, " ≤ ", ub)
            dm[i] = String(take!(cell))
        elseif !isnothing(lb_data) && isnothing(ub_data)
            lb = sep(lb_data[i], 0, 0)
            cell = IOBuffer()
            print(cell, pe, " ≥ ", lb)
            dm[i] = String(take!(cell))
        else
            ubv = ub_data[i]
            lbv = lb_data[i]
            aub = lbv < ubv ? ubv : lbv
            alb = lbv < ubv ? lbv : ubv
            if aub != ubv
                println("[", fd.label, "] This should never print. UB issue ", aub, " ≠ ", ubv)
            elseif alb != lbv
                println("[", fd.label, "] This should never print. LB issue ", alb, " ≠ ", lbv)
            end
            ub = sep(aub, 0, 0)
            lb = sep(alb, 0, 0)
            cell = IOBuffer()
            print(cell, lb, " ≤ ", pe, " ≤ ", ub)
            dm[i] = String(take!(cell))
        end
    end

    title_io = IOBuffer()
    println(title_io, "Forecasted Data: ", fd.label)
    if !isnothing(fd.cl)
        println(title_io, "Prediction Interval Level: ", round(fd.cl * 100, digits=0), "%")
    end
    if !isnothing(fd.uncertainty)
        println(title_io, "Uncertainty Accounted For: ", fd.uncertainty == UM_INNOVATION_ONLY ? "Random Walk Innovation Only" : "Random Walk Innovation + Drift")
    end

    if !isnothing(ub_data) && !isnothing(lb_data)
        println(title_io, "Prediction Interval Format: Lower Bound ≤ Point Forecast ≤ Upper Bound")
    elseif !isnothing(ub_data) && isnothing(lb_data)
        println(title_io, "Prediction Interval Format: Point Forecast ≤ Upper Bound")
    elseif !isnothing(lb_data) && isnothing(ub_data)
        println(title_io, "Prediction Interval Format:  Point Forecast ≥ Lower Bound")
    end
    if rescale
        println(title_io, "Data Units: ", "×10^(", Int(scale_exp), ")")
    end
    title = String(take!(title_io))

    headers = map(Symbol, fd.forecast.periods.values)
    rows = fd.forecast.ages.values
    rt = "Age\\Year"

    fd_config = table_config(io, title=title, rows=rows, row_label_title=rt, headers=headers, alignment=:c, hlines=:all)

    pretty_table(io, dm; fd_config...)

end

function Base.show(io::IO, t::MIME"text/plain", mf::ModelForecasts)
    ds = displaysize(io)
    width = ds[2]

    line = repeat("=", width)

    println(io, line)
    println(io, "Model Forecasts")
    println(io, line)
    println(io)
    show(io, t, mf.kappas)
    println(io, line)
    println(io)
    show(io, t, mf.rates)
    println(io)
    println(io, line)
    println(io)
    show(io, t, mf.deaths)
    println(io)
    println(io, line)
    println(io)
    show(io, t, mf.lifespans)
    println(io)
    println(io, line)


end


center(fd::ForecastedData{ParameterSet{Float64}}) = fd.forecast.values

lower(fd::ForecastedData{ParameterSet{Float64}}) = fd.lower.values

upper(fd::ForecastedData{ParameterSet{Float64}}) = fd.upper.values

center(fd::ForecastedData{AgePeriodData{Float64}}) = fd.forecast.values

lower(fd::ForecastedData{AgePeriodData{Float64}}) = fd.lower.values

upper(fd::ForecastedData{AgePeriodData{Float64}}) = fd.upper.values


Center(fd::ForecastedData{ParameterSet{Float64}}) = fd.forecast

Lower(fd::ForecastedData{ParameterSet{Float64}}) = fd.lower

Upper(fd::ForecastedData{ParameterSet{Float64}}) = fd.upper


kappas(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = [center, lower, upper][Int(fl)](mf.kappas)

rates(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = [center, lower, upper][Int(fl)](mf.rates)

logrates(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = log.([center, lower, upper][Int(fl)](mf.rates))

deaths(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = [center, lower, upper][Int(fl)](mf.deaths)

lifespans(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = [center, lower, upper][Int(fl)](mf.lifespans)



Kappas(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = [Center, Lower, Upper][Int(fl)](mf.kappas)

Rates(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = [Center, Lower, Upper][Int(fl)](mf.rates)

Logrates(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = log.([Center, Lower, Upper][Int(fl)](mf.rates))

Deaths(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = [Center, Lower, Upper][Int(fl)](mf.deaths)

Lifespans(mf::ModelForecasts, fl::ForecastLevel=FL_CENTER) = [Center, Lower, Upper][Int(fl)](mf.lifespans)

