
export maswe

function maswe(model::MortalityModel, forecasts::ModelForecasts; max_lookahead::Int=20, variable::ForecastVariable=FV_LOGRATE)

    exposure_weights = weights(vec(mean(model.exposures.fit.data, dims=2)))

    if variable == FV_LOGRATE
        train = model.logrates.fit.data
        test = model.logrates.test.data
        prediction = log.(forecasts.rates.forecast.data)
    elseif variable == FV_RATE
        train = model.rates.fit.data
        test = model.rates.test.data
        prediction = forecasts.rates.forecast.data
    elseif variable == FV_LE
        train = model.expectedlifetimes.fit.data
        test = model.expectedlifetimes.test.data
        prediction = forecasts.expectedlifetimes.forecast.data
    else
        train = model.deaths.fit.data
        test = model.deaths.test.data
        prediction = model.exposures.test.data .* forecasts.rates.forecast.data
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