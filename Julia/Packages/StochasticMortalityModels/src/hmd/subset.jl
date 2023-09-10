export subset

function subset(df::DataFrame, years::Vector{Int64}, ages::Vector{Int64})
    filter(row -> in(row.Year, years) && in(row.Age, ages), df)
end

function subset(df::DataFrame, years::AbstractRange, ages::AbstractRange)
    allyears = collect(years)
    allages = collect(ages)

    subset(df, allyears, allages)
end


function subset(df::DataFrame; kwargs...)
    arg_keys = keys(kwargs)
    has_min_year = in(:min_year, arg_keys)
    has_max_year = in(:max_year, arg_keys)
    has_min_age = in(:min_age, arg_keys)
    has_max_age = in(:max_age, arg_keys)
    has_year_range = in(:year_range, arg_keys)
    has_age_range = in(:age_range, arg_keys)
    has_years = in(:years, arg_keys)
    has_ages = in(:ages, arg_keys)

    df_max_age = maximum(df.Age)
    df_min_age = minimum(df.Age)
    df_min_year = minimum(df.Year)
    df_max_year = maximum(df.Year)

    age_range = has_age_range ?
                kwargs[:age_range] :
                ((
        has_min_age ?
        kwargs[:min_age] :
        df_min_age
    ):(
        has_max_age ?
        kwargs[:max_age] :
        df_max_age
    ))
    year_range = has_year_range ?
                 kwargs[:year_range] :
                 ((
        has_min_year ?
        kwargs[:min_year] :
        df_min_year
    ):(
        has_max_year ?
        kwargs[:max_year] :
        df_max_year
    ))

    years = has_years ? kwargs[:years] : collect(year_range)
    ages = has_ages ? kwargs[:ages] : collect(age_range)

    subset(df, years, ages)
end
