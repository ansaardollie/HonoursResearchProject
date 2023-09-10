export extract

function extract(df::DataFrame, col::Symbol)

    years = unique(df.Year)
    ages = unique(df.Age)
    data = df[!, col]

    dm = Matrix{Float64}(reshape(data, length(ages), length(years)))

    return dm
end

