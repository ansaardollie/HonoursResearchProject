
export maswe, residuals, forecast_errors

function maswe(model::MortalityModel, forecasts::ModelForecasts; max_lookahead::Int=20, variable::ForecastVariable=FV_LOGRATE)
    n_ages = length(ages(model))
    n_years_test = length(years(model, DS_TEST))
    age_weights = variable == FV_LEAB ? weights([1]) : weights(vec(mean(exposures(model, DS_TRAIN), dims=2)))

    if variable == FV_LEAB
        func = Lifespans
    else
        func = [LogRates, Rates, Lifespans, Deaths][Int(variable)]
    end

    first_year = years(model, DS_COMPLETE)[begin]
    last_train_year = years(model, DS_TRAIN)[end]
    train = func(model, DS_COMPLETE)[:, first_year:last_train_year].values
    test = func(model, DS_TEST).values

    if variable == FV_LEAB
        train = reshape(train[1, :], 1, size(train, 2))
        test = reshape(test[1, :], 1, size(test, 2))
    end

    if variable == FV_LOGRATE
        prediction = logrates(forecasts)
    elseif variable == FV_RATE
        prediction = rates(forecasts)
    elseif variable == FV_MRL
        prediction = lifespans(forecasts)
    elseif variable == FV_DEATHS
        prediction = deaths(forecasts)
    else
        prediction = reshape(lifespans(forecasts)[1, :], 1, n_years_test)
    end

    results = Vector{Float64}(undef, max_lookahead)
    for l in 1:max_lookahead
        naive_true = train[:, (begin+l):end]
        naive_forecast = train[:, begin:(end-l)]
        naive_mvae = abs.(naive_true - naive_forecast)
        naive_waae = mean(naive_mvae, age_weights, dims=1)

        fc_errors = abs.(test - prediction)
        fc_waae = mean(fc_errors, age_weights, dims=1)

        mae_naive = mean(naive_waae)
        mae_fc = mean(fc_waae)
        results[l] = mae_fc / mae_naive
    end

    return results

end


function residuals(model::MortalityModel; variable::ForecastVariable=FV_LOGRATE)
    if variable == FV_LEAB
        throw(ArgumentError("Can't use Life Expectancy At Birth"))
    end
    func = [logrates, rates, lifespans, deaths][Int(variable)]

    actual = func(model, DS_TRAIN)

    if variable == FV_LOGRATE
        fitted = mxt_hat(model.parameters, log_scale=true).values
    elseif variable == FV_RATE
        fitted = mxt_hat(model.parameters, log_scale=false).values
    elseif variable == FV_MRL
        mxt = mxt_hat(model.parameters, log_scale=false).values
        x = ages(model, DS_TRAIN)
        t = years(model, DS_TRAIN)
        fitted = lexpectancies(mxt, x, t; gender=sex(model), at_age=x)
    else
        fitted = fitted_deaths(model).values
    end

    residuals = actual .- fitted

    return residuals
end

function forecast_errors(model::MortalityModel, forecasts::ModelForecasts; variable::ForecastVariable=FV_LOGRATE)
    if variable == FV_LEAB
        throw(ArgumentError("Can't use Life Expectancy At Birth"))
    end
    func = [logrates, rates, lifespans, deaths][Int(variable)]

    actual = func(model, DS_TEST)

    if variable == FV_LOGRATE
        # fitted = mxt_hat(model.parameters, log_scale=true).values
        fcv = logrates(forecasts)
    elseif variable == FV_RATE
        fcv = rates(forecasts)
    elseif variable == FV_MRL
        fcv = lifespans(forecasts)
    else
        fcv = deaths(forecasts)
    end

    errors = actual .- fcv

    # return AgePeriodData{Float64}(MDC_OTHER, MDS_CALCULATED, "Forecast Error $variable", errors, Ages(model), Years(model), 6, false, Nothing)
    return errors
end