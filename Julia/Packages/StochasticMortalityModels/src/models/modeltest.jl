
export maswe

function maswe(model::MortalityModel, forecasts::ModelForecasts; max_lookahead::Int=20, variable::ForecastVariable=FV_LOGRATE)

    exposure_weights = weights(vec(mean(model.exposures.fit.values, dims=2)))

    if variable == FV_LOGRATE
        train = model.logrates.fit.values
        test = model.logrates.test.values
        prediction = log.(forecasts.rates.forecast.values)
    elseif variable == FV_RATE
        train = model.rates.fit.values
        test = model.rates.test.values
        prediction = forecasts.rates.forecast.values
    elseif variable == FV_LE
        train = model.expectedlifetimes.fit.values
        test = model.expectedlifetimes.test.values
        prediction = forecasts.expectedlifetimes.forecast.values
    else
        train = model.deaths.fit.values
        test = model.deaths.test.values
        prediction = model.exposures.test.values .* forecasts.rates.forecast.values
    end

    results = Vector{Float64}(undef, max_lookahead)
    for l in 1:max_lookahead
        naive_true = train[:, (begin+l):end]
        naive_forecast = train[:, begin:(end-l)]
        naive_mvae = abs.(naive_true - naive_forecast)
        naive_waae = mean(naive_mvae, exposure_weights, dims=1)

        fc_errors = abs.(test - prediction)
        fc_waae = mean(fc_errors, exposure_weights, dims=1)

        mae_naive = mean(naive_waae)
        mae_fc = mean(fc_waae)
        results[l] = mae_fc / mae_naive
    end

    return results

end