include("models.jl")
function chi_squared_test(model::MortalityModel, forecasts::ModelForecasts)

    observed = model.deaths.test.data
    expected = model.exposures.test.data .* forecasts.rates.forecast.data

    st_devs = (observed .- expected) / sqrt.(expected)
end