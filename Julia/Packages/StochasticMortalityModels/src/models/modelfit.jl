export fitrange!,
    basefit!,
    LinearAlgebra,
    adjustkappa_lc!,
    adjustkappa_lc_demography!,
    adjustkappa_lc_demography2!,
    fitted_deaths,
    adjustkappa_lm!,
    adjustkappa_bms!,
    fit!,
    choose_period!


function fitrange!(model::MortalityModel; years::SecondaryRangeSelection=nothing, ages::SecondaryRangeSelection=nothing)

    if isnothing(years) && isnothing(ages)
        throw(ArgumentError("Must provide range to fit for"))
    end

    yr = @match years begin
        ::Nothing => model.ranges.all.years
        ::DataRange => years
        ::AbstractArray{Int} => yearrange(years)
        _ => throw(ArgumentError("Cant pass $(typeof(years)), must be `nothing` a `DataRange object or `Vector{Int}`"))
    end

    ar = @match ages begin
        ::Nothing => model.ranges.all.ages
        ::DataRange => ages
        ::AbstractArray{Int} => agerange(ages)
        _ => throw(ArgumentError("Cant pass $(typeof(years)), must be `nothing` a `DataRange object or `Vector{Int}`"))
    end

    exposures = model.exposures
    @reset exposures.fit = model.exposures.all[ar, yr]
    deaths = model.deaths
    @reset deaths.fit = model.deaths.all[ar, yr]
    approximatedeaths = model.approximatedeaths
    @reset approximatedeaths.fit = model.approximatedeaths.all[ar, yr]
    rates = model.rates
    @reset rates.fit = model.rates.all[ar, yr]
    logrates = model.logrates
    @reset logrates.fit = model.logrates.all[ar, yr]
    expectedlifetimes = model.expectedlifetimes
    @reset expectedlifetimes.fit = model.expectedlifetimes.all[ar, yr]
    ranges::Stratified{AgeYearRange} = Stratified(
        model.ranges.all,
        (years=yr, ages=ar)::AgeYearRange,
        model.ranges.test
    )


    model.exposures = exposures
    model.deaths = deaths
    model.approximatedeaths = approximatedeaths
    model.rates = rates
    model.logrates = logrates
    model.expectedlifetimes = expectedlifetimes
    model.ranges = ranges

    newparams = ModelParameters((years=yr, ages=ar))

    model.parameters = newparams
    return model
end


function fitted_deaths(model::MortalityModel)
    params = model.parameters
    exposures = model.exposures.fit
    return fitted_deaths(params, exposures)
end

function basefit!(model::MortalityModel; constrain::Bool=true)
    α = mapslices(mean, model.logrates.fit.data, dims=2)
    alpha = ParameterSet("α(x)", vec(α), model.ranges.fit.ages)
    Zxt = model.logrates.fit.data .- α

    Zxt_svd = svd(Zxt)
    β = Zxt_svd.U[:, 1]
    beta = ParameterSet("β(x)", vec(β), model.ranges.fit.ages)

    κ = Zxt_svd.V[:, 1] * Zxt_svd.S[1]
    kappa = ParameterSet("κ(t)", vec(κ), model.ranges.fit.years)

    mp = ModelParameters(alpha, beta, kappa)
    if constrain
        constrain!(mp; mode=model.calculationmode)
    end

    model.parameters = mp
end

function basefit!(logrates::Matrix{Float64}; constrain_julia::Bool=true, constrain_demography::Bool=false)
    α = vec(mapslices(mean, logrates, dims=2))

    Zxt = logrates .- α

    Zxt_svd = svd(Zxt)
    β = vec(Zxt_svd.U[:, 1])


    κ = vec(Zxt_svd.V[:, 1] * Zxt_svd.S[1])
    if constrain_julia
        c1 = sum(β)
        c2 = mean(κ)
        α = α .+ (c2 .* β)
        β = β ./ c1
        κ = c1 .* (κ .- c2)
    elseif constrain_demography
        c1 = sum(β)
        β = β ./ c1
        κ = c1 .* κ
    end

    return (alpha=α, beta=β, kappa=κ)
end

function adjustkappa_lc!(model::MortalityModel)
    parameters = deepcopy(model.parameters)

    deaths = model.calculationmode == CM_JULIA ?
             model.deaths.fit.data :
             model.approximatedeaths.fit.data

    betas = parameters.betas.values
    solve_func(fv, jv, kappas) = begin
        @reset parameters.kappas.values = kappas
        fdxt = fitted_deaths(parameters, model.exposures.fit).data
        jdxt = fdxt .* betas
        func_vals = vec(mapslices(sum, fdxt, dims=1)) - vec(mapslices(sum, deaths, dims=1))
        jacob_diag = vec(mapslices(sum, jdxt, dims=1))

        if !isnothing(fv)
            copyto!(fv, func_vals)
        end

        if !isnothing(jv)
            k = length(jacob_diag)
            for i in 1:k, j in 1:k
                jv[i, j] = (i == j) ? jacob_diag[i] : 0
            end
        end
    end

    solve_out = nlsolve(only_fj!(solve_func), parameters.kappas.values)
    @reset parameters.kappas.values = solve_out.zero

    if model.calculationmode == CM_JULIA
        model.parameters = constrain_julia!(parameters)
    else
        model.parameters = parameters
    end


end

function adjustkappa_lm!(model::MortalityModel)
    M = length(model.ranges.fit.years)

    parameters = deepcopy(model.parameters)
    ages = model.ranges.fit.ages.values
    years = model.ranges.fit.years.values
    start_age = ages[1]
    obs_e0 = fill(0, M)

    if model.calculationmode == CM_JULIA
        obs_e0 = vec(model.expectedlifetimes.fit[start_age, :].data)
    else
        mxt = model.rates.fit.data
        el = mapslices(mx -> expected_lifetime(mx, ages; sex=model.population.sex, at_age=[start_age], mode=model.calculationmode), mxt, dims=1)
        obs_e0 = vec(el)
    end

    obs_e0_ps = ParameterSet("Expected Lifetimes", obs_e0, model.ranges.fit.years)

    opt_kappas = Vector{Float64}()
    for year in years
        kappas = parameters.kappas[year].values
        alphas = parameters.alphas.values
        betas = parameters.betas.values

        e0 = obs_e0_ps[year].values[1]
        solve_func!(F, kv) = begin
            kappa = kv[1]
            mxt = exp.(alphas + (betas * kappa))
            le = expected_lifetime(mxt, ages; sex=model.population.sex, mode=model.calculationmode, at_age=[start_age])
            F[1] = le[1] - e0
        end

        solve_out = nlsolve(solve_func!, kappas, autodiff=:forward, method=:newton)

        push!(opt_kappas, solve_out.zero[1])
    end

    @reset parameters.kappas.values = opt_kappas

    model.parameters = parameters
end

function deviance(obs, fit)

    if obs == 0 || fit == 0
        return 0
    else
        return 2 * (obs * log(obs / fit) - (obs - fit))
    end
end


function adjustkappa_bms!(model::MortalityModel; constrain::Bool=true)
    parameters = deepcopy(model.parameters)
    years = model.ranges.fit.years.values
    deaths = model.calculationmode == CM_JULIA ? model.deaths.fit.data : model.approximatedeaths.fit.data
    exposures = model.exposures.fit.data

    alphas = parameters.alphas.values
    betas = parameters.betas.values
    init_kappas = parameters.kappas.values

    opt_kappas = Vector{Float64}()
    for i in eachindex(years)
        year = years[i]
        dt_byx = vec(deaths[:, i])
        et_byx = vec(exposures[:, i])
        k0 = init_kappas[i]

        obj_func(kv) = begin
            kappa = kv[1]
            mx = exp.(alphas + betas * kappa)
            fdx = et_byx .* mx
            dev = deviance.(dt_byx, fdx)
            return sum(dev)
        end

        grad_func(G, kv) = begin
            kappa = kv[1]
            mx = exp.(alphas + betas * kappa)
            fdx = et_byx .* mx
            gterms = 2 * (betas .* (fdx - dt_byx))
            G[1] = sum(gterms)
        end

        opt_result = optimize(obj_func, grad_func, [k0])
        push!(opt_kappas, opt_result.minimizer[1])
    end


    @reset parameters.kappas.values = opt_kappas

    if constrain
        parameters = constrain!(parameters; mode=model.calculationmode)
    end

    model.parameters = parameters
end


function adjustkappa_bms!(deaths::Matrix{Float64}, exposures::Matrix{Float64}, alphas::Vector{Float64}, betas::Vector{Float64}, kappas::Vector{Float64}; constrain_julia::Bool=true)

    init_kappas = kappas

    opt_kappas = Vector{Float64}()
    for i in eachindex(init_kappas)

        dt_byx = vec(deaths[:, i])
        et_byx = vec(exposures[:, i])
        k0 = init_kappas[i]

        obj_func(kv) = begin
            kappa = kv[1]
            mx = exp.(alphas + betas * kappa)
            fdx = et_byx .* mx
            dev = deviance.(dt_byx, fdx)
            return sum(dev)
        end

        grad_func(G, kv) = begin
            kappa = kv[1]
            mx = exp.(alphas + betas * kappa)
            fdx = et_byx .* mx
            gterms = 2 * (betas .* (fdx - dt_byx))
            G[1] = sum(gterms)
        end

        opt_result = optimize(obj_func, grad_func, [k0])
        push!(opt_kappas, opt_result.minimizer[1])
    end

    α = alphas
    β = betas
    κ = opt_kappas

    if constrain_julia
        c1 = sum(β)
        c2 = mean(κ)
        α = α .+ (c2 .* β)
        β = β ./ c1
        κ = c1 .* (κ .- c2)
    end

    return (alpha=α, beta=β, kappa=κ)
end



function mt_adjustkappa_lm!(model::MortalityModel)
    M = length(model.ranges.fit.years)

    parameters = deepcopy(model.parameters)
    ages = model.ranges.fit.ages.values
    years = model.ranges.fit.years.values
    start_age = ages[1]
    obs_e0 = fill(0, M)

    if model.calculationmode == CM_JULIA
        obs_e0 = vec(model.expectedlifetimes.fit[start_age, :].data)
    else
        mxt = model.rates.fit.data
        el = mapslices(mx -> expected_lifetime(mx, ages; sex=model.population.sex, at_age=[start_age], mode=model.calculationmode), mxt, dims=1)
        obs_e0 = vec(el)
    end

    obs_e0_ps = ParameterSet("Expected Lifetimes", obs_e0, model.ranges.fit.years)

    opt_kappas = Dict{Int,Float64}()
    res_channel = Channel(c -> begin
        @threads for year in years
            kappas = parameters.kappas[year].values
            alphas = parameters.alphas.values
            betas = parameters.betas.values

            e0 = obs_e0_ps[year].values[1]
            solve_func!(F, kv) = begin
                kappa = kv[1]
                mxt = exp.(alphas + (betas * kappa))
                le = expected_lifetime(mxt, ages; sex=model.population.sex, mode=model.calculationmode, at_age=[start_age])
                F[1] = le[1] - e0
            end

            solve_out = nlsolve(solve_func!, kappas, autodiff=:forward, method=:newton)

            put!(c, (year=year, kappa=solve_out.zero[1]))
        end
    end)


    for res in res_channel
        year = res.year
        kappa = res.kappa
        opt_kappas[year] = kappa
    end

    @reset parameters.kappas.values = [opt_kappas[k] for k in sort(collect(keys(opt_kappas)))]

    model.parameters = parameters
end


function mt_adjustkappa_bms!(model::MortalityModel)
    parameters = deepcopy(model.parameters)
    years = model.ranges.fit.years.values
    deaths = model.calculationmode == CM_JULIA ? model.deaths.fit.data : model.approximatedeaths.fit.data
    exposures = model.exposures.fit.data

    alphas = parameters.alphas.values
    betas = parameters.betas.values
    init_kappas = parameters.kappas.values

    opt_kappas = Dict{Int,Float64}()
    res_channel = Channel(c -> begin
        @threads for i in eachindex(years)
            year = years[i]
            dt_byx = vec(deaths[:, i])
            et_byx = vec(exposures[:, i])
            k0 = init_kappas[i]

            obj_func(kv) = begin
                kappa = kv[1]
                mx = exp.(alphas + betas * kappa)
                fdx = et_byx .* mx
                dev = deviance.(dt_byx, fdx)
                return sum(dev)
            end

            grad_func(G, kv) = begin
                kappa = kv[1]
                mx = exp.(alphas + betas * kappa)
                fdx = et_byx .* mx
                gterms = 2 * (betas .* (fdx - dt_byx))
                G[1] = sum(gterms)
            end

            opt_result = optimize(obj_func, grad_func, [k0])
            put!(c, (year=year, kappa=opt_result.minimizer[1]))
        end
    end)

    for res in res_channel
        year = res.year
        kappa = res.kappa
        opt_kappas[year] = kappa
    end



    @reset parameters.kappas.values = [opt_kappas[k] for k in sort(collect(keys(opt_kappas)))]

    if model.calculationmode == CM_JULIA
        parameters = constrain_julia!(parameters)
    end

    model.parameters = parameters
end


function fit!(model::MortalityModel; constrain::Bool=true, multithreaded::Bool=false, choose_period::Bool=false)

    if choose_period
        choose_period!(model)
    end

    basefit!(model, constrain=constrain)

    if model.variant.adjustment == AC_DXT
        multithreaded ? mt_adjustkappa_bms!(model) : adjustkappa_bms!(model)
    elseif model.variant.adjustment == AC_E0
        multithreaded ? mt_adjustkappa_lm!(model) : adjustkappa_lm!(model)
    else
        adjustkappa_lc!(model)
    end

    return model

end


function calculate_deviance_statistic(deaths::Matrix{Float64}, exposures::Matrix{Float64}, alphas::Vector{Float64}, betas::Vector{Float64}, kappas::Vector{Float64})
    X = size(deaths, 1)
    T = size(deaths, 2)
    mxt = exp.(reshape(alphas, X, 1) .+ reshape(betas, X, 1) * reshape(kappas, 1, T))
    fdxt = exposures .* mxt
    dev_xt = deviance.(deaths, fdxt)

    mk = mean(kappas)
    slope = (kappas[end] - kappas[begin]) / (T - 1)
    mt = 1 + (T / 2)

    t = collect(1:T)
    lf_kappas = mk .+ slope .* (t .- mt)

    lf_mxt = exp.(reshape(alphas, X, 1) .+ reshape(betas, X, 1) * reshape(lf_kappas, 1, T))
    lf_fdxt = exposures .* lf_mxt
    dev_lfxt = deviance.(deaths, lf_fdxt)

    dev_base = sum(dev_xt)
    dev_total = sum(dev_lfxt)

    df_base = ((X - 1) * (T - 2))
    df_total = (X * (T - 2))

    R_s = (dev_total / df_total) / (dev_base / df_base)

    return R_s

end


function choose_period!(model::MortalityModel)
    all_years = model.ranges.all.years.values
    fit_years = model.ranges.fit.years.values
    test_years = model.ranges.test.years.values

    if length(all_years) != length(test_years) || !all(all_years .== test_years)
        years = collect(all_years[begin]:(test_years[begin]-1))
    else
        years = model.ranges.fit.years.values
    end


    m = years[end]
    potential_starts = years[1:end-2]
    Sl = length(potential_starts)

    cj = model.calculationmode == CM_JULIA
    cd = model.calculationmode == CM_DEMOGRAPHY

    R_statistics = Matrix{Float64}(undef, Sl, 2)

    for i in eachindex(potential_starts)
        S = potential_starts[i]
        lmxt = model.logrates.fit[:, S:m].data
        dxt = model.deaths.fit[:, S:m].data
        Ext = model.exposures.fit[:, S:m].data
        (alpha, beta, kappa) = basefit!(lmxt; constrain_julia=cj, constrain_demography=cd)

        (alpha, beta, kappa) = adjustkappa_bms!(dxt, Ext, alpha, beta, kappa; constrain_julia=false)

        R_S = calculate_deviance_statistic(dxt, Ext, alpha, beta, kappa)

        R_statistics[i, :] = [S, R_S]
    end

    ordered_idx = sortperm(R_statistics[:, 2])
    start_year = Int(R_statistics[ordered_idx[1], 1])
    yr = yearrange(start_year:m)

    fitrange!(model; years=yr)

end


function variation_explained(model::MortalityModel)
    logrates = model.logrates.all.data
    α = mapslices(mean, logrates, dims=2)
    centered = logrates .- α

    σ = svd(centered).S

    ∑σ² = sum(σ .^ 2)

    pve = cumsum(σ .^ 2) ./ ∑σ²

    return pve
end

