
export maswe

function maswe(model::MortalityModel, forecasts::ModelForecasts; max_lookahead::Int=20, variable::ForecastVariable=FV_LOGRATE)

    exposure_weights = weights(vec(mean(exposures(model, DS_TRAIN), dims=2)))

    func = [logrates, rates, lifespans, deaths][Int(variable)]

    train = func(model, DS_TRAIN)
    test = func(model, DS_TEST)

    if variable == FV_LOGRATE
        prediction = logrates(forecasts)
    elseif variable == FV_RATE
        prediction = rates(forecasts)
    elseif variable == FV_LE
        prediction = lifespans(forecasts)
    else
        prediction = deaths(forecasts)
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