export fitrange!,
    basefit!,
    fitted_deaths,
    fit!,
    choose_period!,
    deviance,
    lc_adjust,
    lc_adjust!,
    lm_adjust,
    lm_adjust!,
    bms_adjust,
    bms_adjust!,
    variation_explained,
    gen_kopt_bms,
    gen_kopt_lc,
    gen_kopt_lm

SecondaryRangeSelection = Union{Nothing,DataRange,AbstractVector{Int}}

function fitrange!(model::MortalityModel; yvals::SecondaryRangeSelection=nothing, avals::SecondaryRangeSelection=nothing)

    if isnothing(yvals) && isnothing(avals)
        throw(ArgumentError("Must provide range to fit for"))
    end

    yr::DataRange = @match yvals begin
        ::Nothing => Years(model, DS_COMPLETE)
        ::DataRange => yvals
        ::AbstractArray{Int} => yearrange(yvals)
        _ => throw(ArgumentError("Cant pass $(typeof(yvals)), must be `nothing` a `DataRange object or `Vector{Int}`"))
    end

    ar::DataRange = @match avals begin
        ::Nothing => Ages(model, DS_COMPLETE)
        ::DataRange => avals
        ::AbstractArray{Int} => agerange(avals)
        _ => throw(ArgumentError("Cant pass $(typeof(yvals)), must be `nothing` a `DataRange object or `Vector{Int}`"))
    end


    yidx = collect(map(year -> Years(model, DS_COMPLETE).indexer[year], yr.values))
    aidx = collect(map(age -> Ages(model, DS_COMPLETE).indexer[age], ar.values))
    train::MortalityData = model.data.train

    @reset train.exposures = exposures(model, DS_COMPLETE)[aidx, yidx]
    @reset train.deaths = deaths(model, DS_COMPLETE)[aidx, yidx]
    @reset train.rates = rates(model, DS_COMPLETE)[aidx, yidx]
    @reset train.logrates = logrates(model, DS_COMPLETE)[aidx, yidx]
    @reset train.lifespans = lifespans(model, DS_COMPLETE)[aidx, yidx]

    @reset train.ages = ar.values
    @reset train.years = yr.values

    model.data = Stratified{MortalityData}(
        model.data.complete,
        train,
        model.data.test,
        model.data.label
    )

    model.ranges = Stratified{AgeYearRange}(
        model.ranges.complete,
        (years=yr, ages=ar),
        model.ranges.test,
        model.ranges.label
    )

    newparams = ModelParameters((years=yr, ages=ar))

    model.parameters = newparams
    return model
end


function fitted_deaths(model::MortalityModel, version::ParameterVersion=PV_ADJUSTED)
    params = version == PV_ADJUSTED ? model.parameters : model.unadjusted_parameters
    return fitted_deaths(params, Exposures(model))
end

function basefit!(model::MortalityModel; constrain::Bool=true)
    ar = Ages(model)
    yr = Years(model)
    lr = logrates(model)
    α = mapslices(lmt -> mean(lmt[.!isinf.(lmt)]), lr, dims=2)
    alpha = ParameterSet("α(x)", vec(α), ar)
    Zxt = lr .- α

    Zxt_svd = svd(Zxt)
    β = Zxt_svd.U[:, 1]
    beta = ParameterSet("β(x)", vec(β), ar)

    κ = Zxt_svd.V[:, 1] * Zxt_svd.S[1]
    kappa = ParameterSet("κ(t)", vec(κ), yr)

    mp = ModelParameters(alpha, beta, kappa)
    ump = ModelParameters(alpha, beta, kappa)
    if constrain
        constrain!(mp; mode=calculation(model))
        constrain!(ump; mode=calculation(model))
    end

    model.parameters = mp
    model.unadjusted_parameters = ump
end

function basefit!(logrates::Matrix{Float64}; constrain::Bool=true, cmode::CalculationChoice=CC_JULIA)
    α = vec(mapslices(lmt -> mean(lmt[.!isinf.(lmt)]), logrates, dims=2))

    Zxt = logrates .- α

    Zxt_svd = svd(Zxt)
    β = vec(Zxt_svd.U[:, 1])


    κ = vec(Zxt_svd.V[:, 1] * Zxt_svd.S[1])
    if constrain && cmode == CC_JULIA
        c1 = sum(β)
        c2 = mean(κ)
        α = α .+ (c2 .* β)
        β = β ./ c1
        κ = c1 .* (κ .- c2)
    elseif constrain && cmode == CC_DEMOGRAPHY
        c1 = sum(β)
        β = β ./ c1
        κ = c1 .* κ
    end

    return (alphas=α, betas=β, kappas=κ)
end

function deviance(obs, fit)

    if obs == 0 || fit == 0
        return 0
    else
        return 2 * (obs * log(obs / fit) - (obs - fit))
    end
end

AdjustmentInfoLC = NamedTuple{(:init, :deaths, :exposures),Tuple{Float64,Vector{Float64},Vector{Float64}}}

AdjustmentInfoLM = NamedTuple{(:init, :e0),Tuple{Float64,Float64}}

AdjustmentInfoBMS = NamedTuple{(:init, :deaths, :exposures),Tuple{Float64,Vector{Float64},Vector{Float64}}}

function gen_kopt_lc(
    ax::Vector{Float64},
    bx::Vector{Float64}
)

    n = length(ax)
    gen_func(info::AdjustmentInfoLC) = begin
        mx_f = Vector{Float64}(fill(0, n))
        mx_j = Vector{Float64}(fill(0, n))
        fdx_f = Vector{Float64}(fill(0, n))
        fdx_j = Vector{Float64}(fill(0, n))

        solve_func!(F, kv) = begin
            kappa = kv[1]
            copyto!(mx_f, exp.(ax + bx * kappa))
            copyto!(fdx_f, info.exposures .* mx_f)

            F[1] = sum(fdx_f .- info.deaths)
        end

        jacob_func!(G, kv) = begin
            kappa = kv[1]
            copyto!(mx_j, exp.(ax + bx * kappa))
            copyto!(fdx_j, info.exposures .* mx_j)
            G[1] = sum(fdx_j .* bx)
        end

        solve_result = nlsolve(solve_func!, jacob_func!, [info.init])

        return solve_result.zero[1]
    end

    return gen_func
end

function gen_kopt_lm(
    ax::Vector{Float64},
    bx::Vector{Float64},
    x::Vector{Int},
    gender::Sex,
    start_age::Int,
    mode::CalculationChoice=CC_JULIA
)

    n = length(ax)
    gen_func(info::AdjustmentInfoLM) = begin

        solve_func!(F, kv) = begin
            kappa = kv[1]
            mx = exp.(ax + bx * kappa)
            fe0 = expected_lifetime(mx, x; sex=gender, at_age=[start_age], mode=mode)
            F[1] = fe0[1] - info.e0
        end


        solve_result = nlsolve(solve_func!, [info.init], autodiff=:forward)

        return solve_result.zero[1]
    end

    return gen_func
end

function gen_kopt_bms(
    ax::Vector{Float64},
    bx::Vector{Float64}
)

    n = length(ax)
    gen_func(info::AdjustmentInfoBMS) = begin
        mx_f = Vector{Float64}(fill(0, n))
        mx_g = Vector{Float64}(fill(0, n))
        fdx_f = Vector{Float64}(fill(0, n))
        fdx_g = Vector{Float64}(fill(0, n))
        dev_f = Vector{Float64}(fill(0, n))
        grad_t = Vector{Float64}(fill(0, n))
        obj_func(kv) = begin
            kappa = kv[1]
            copyto!(mx_f, exp.(ax + bx * kappa))
            copyto!(fdx_f, info.exposures .* mx_f)
            copyto!(dev_f, deviance.(info.deaths, fdx_f))
            return sum(dev_f)
        end

        grad_func(G, kv) = begin
            kappa = kv[1]
            copyto!(mx_g, exp.(ax + bx * kappa))
            copyto!(fdx_g, info.exposures .* mx_g)
            copyto!(grad_t, 2 * (bx .* (fdx_g - info.deaths)))
            G[1] = sum(grad_t)
        end

        opt_result = optimize(obj_func, grad_func, [info.init])

        return opt_result.minimizer[1]
    end

    return gen_func
end

function lc_adjust(
    ax::Vector{Float64},
    bx::Vector{Float64},
    kt::Vector{Float64},
    t::Vector{Int},
    dxt::Matrix{Float64},
    Ext::Matrix{Float64};
    constrain::Bool=true,
    cmode::CalculationChoice=CC_JULIA
)
    data = Dict{Int,AdjustmentInfoLC}()

    for i in eachindex(t)
        data[t[i]] = (
            init=kt[i],
            deaths=dxt[:, i],
            exposures=Ext[:, i]
        )
    end

    opt_f = gen_kopt_lc(ax, bx)

    output = ThreadsX.map(data) do pair
        return pair[1] => opt_f(pair[2])
    end

    opt_kappas = map(x -> x[2], sort!(output, by=p -> p[1]))

    if constrain && cmode == CC_JULIA
        c1 = sum(bx)
        c2 = mean(opt_kappas)

        α = ax .+ c2 * bx
        β = bx ./ c1
        κ = c1 * (opt_kappas .- c2)
    elseif constrain && cmode == CC_DEMOGRAPHY
        c = sum(bx)
        α = ax
        β = bx ./ c
        κ = c .* opt_kappas
    else
        α = ax
        β = bx
        κ = opt_kappas
    end

    return (alphas=α, betas=β, kappas=κ)
end

function lm_adjust(
    ax::Vector{Float64},
    bx::Vector{Float64},
    kt::Vector{Float64},
    t::Vector{Int},
    x::Vector{Int},
    mxt::Matrix{Float64},
    ext::Matrix{Float64},
    gender::Sex;
    cmode::CalculationChoice=CC_JULIA,
    constrain::Bool=true
)
    data = Dict{Int,AdjustmentInfoLM}()

    start_age = x[1]
    if cmode == CC_JULIA
        obs_e0 = ext[1, :]
    else
        el = mapslices(mx -> expected_lifetime(mx, x; sex=gender, at_age=[start_age], mode=cmode), mxt, dims=1)
        obs_e0 = vec(el)
    end

    for i in eachindex(t)
        data[t[i]] = (
            init=kt[i],
            e0=obs_e0[i])
    end

    solve_f = gen_kopt_lm(ax, bx, x, gender, start_age, cmode)

    output = ThreadsX.map(data) do pair
        return pair[1] => solve_f(pair[2])
    end

    opt_kappas = map(x -> x[2], sort!(output, by=p -> p[1]))

    if constrain && cmode == CC_JULIA
        c1 = sum(bx)
        c2 = mean(opt_kappas)

        α = ax .+ c2 * bx
        β = bx ./ c1
        κ = c1 * (opt_kappas .- c2)
    elseif constrain && cmode == CC_DEMOGRAPHY
        c = sum(bx)
        α = ax
        β = bx ./ c
        κ = c .* opt_kappas
    else
        α = ax
        β = bx
        κ = opt_kappas
    end

    return (alphas=α, betas=β, kappas=κ)
end

function bms_adjust(
    ax::Vector{Float64},
    bx::Vector{Float64},
    kt::Vector{Float64},
    t::Vector{Int},
    dxt::Matrix{Float64},
    Ext::Matrix{Float64};
    constrain::Bool=true,
    cmode::CalculationChoice=CC_JULIA
)
    data = Dict{Int,AdjustmentInfoBMS}()

    for i in eachindex(t)
        data[t[i]] = (
            init=kt[i],
            deaths=dxt[:, i],
            exposures=Ext[:, i]
        )
    end

    opt_f = gen_kopt_bms(ax, bx)

    output = ThreadsX.map(data) do pair
        return pair[1] => opt_f(pair[2])
    end

    opt_kappas = map(x -> x[2], sort!(output, by=p -> p[1]))

    if constrain && cmode == CC_JULIA
        c1 = sum(bx)
        c2 = mean(opt_kappas)

        α = ax .+ c2 * bx
        β = bx ./ c1
        κ = c1 * (opt_kappas .- c2)
    elseif constrain && cmode == CC_DEMOGRAPHY
        c = sum(bx)
        α = ax
        β = bx ./ c
        κ = c .* opt_kappas
    else
        α = ax
        β = bx
        κ = opt_kappas
    end

    return (alphas=α, betas=β, kappas=κ)
end

function lc_adjust!(model::MortalityModel; constrain::Bool=true)
    α = alphas(model)
    β = betas(model)
    κ = kappas(model)
    t = years(model)
    x = ages(model)
    dxt = deaths(model)
    Ext = exposures(model)
    cm = calculation(model)
    (α′, β′, κ′) = lc_adjust(α, β, κ, t, dxt, Ext; constrain=constrain, cmode=cm)
    θ = ModelParameters(α′, β′, κ′, x, t)

    model.parameters = θ
    return model
end

function lm_adjust!(model::MortalityModel; constrain::Bool=true)
    α = alphas(model)
    β = betas(model)
    κ = kappas(model)
    t = years(model)
    x = ages(model)
    mxt = rates(model)
    ext = lifespans(model)
    gender = sex(model)
    cmode = calculation(model)
    (α′, β′, κ′) = lm_adjust(α, β, κ, t, x, mxt, ext, gender; cmode=cmode, constrain=constrain)
    θ = ModelParameters(α′, β′, κ′, x, t)
    model.parameters = θ
    return model
end

function bms_adjust!(model::MortalityModel; constrain::Bool=true)
    α = alphas(model)
    β = betas(model)
    κ = kappas(model)
    t = years(model)
    x = ages(model)
    dxt = deaths(model)
    Ext = exposures(model)
    cm = calculation(model)
    (α′, β′, κ′) = bms_adjust(α, β, κ, t, dxt, Ext; constrain=constrain, cmode=cm)
    θ = ModelParameters(α′, β′, κ′, x, t)
    model.parameters = θ
    return model
end

function calculate_deviance_statistic(
    deaths::Matrix{Float64},
    exposures::Matrix{Float64},
    alphas::Vector{Float64},
    betas::Vector{Float64},
    kappas::Vector{Float64}
)
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

function choose_period!(model::MortalityModel; constrain::Bool=true)
    all_years = years(model, DS_COMPLETE)
    train_years = years(model, DS_TRAIN)
    test_years = years(model, DS_TEST)

    if length(all_years) != length(test_years) || !all(all_years .== test_years)
        yvals = collect(all_years[begin]:(test_years[begin]-1))
    else
        yvals = train_years
    end

    yr = yearrange(yvals)
    m = yvals[end]
    m_idx = yr.indexer[m]
    potential_starts = yvals[1:end-2]
    ps_idx = collect(map(y -> yr.indexer[y], potential_starts))
    Sl = length(potential_starts)

    cm = calculation(model)

    R_statistics = Matrix{Float64}(undef, Sl, 2)

    for i in eachindex(ps_idx)
        S = potential_starts[i]
        S_idx = ps_idx[i]
        lmxt = logrates(model)[:, S_idx:m_idx]
        dxt = deaths(model)[:, S_idx:m_idx]
        Ext = exposures(model)[:, S_idx:m_idx]
        (α, β, κ) = basefit!(lmxt; constrain=constrain, cmode=cm)

        (α, β, κ) = bms_adjust(α, β, κ, collect(S:m), dxt, Ext; constrain=constrain && cm == CC_JULIA, cmode=cm)

        R_S = calculate_deviance_statistic(dxt, Ext, α, β, κ)

        R_statistics[i, :] = [S, R_S]
    end

    ordered_idx = sortperm(R_statistics[:, 2])
    start_year = Int(R_statistics[ordered_idx[1], 1])
    yr = yearrange(start_year:m)

    fitrange!(model; yvals=yr)

    return R_statistics

end

function fit!(model::MortalityModel; constrain::Bool=true, choose_period::Bool=false)

    if choose_period
        choose_period!(model)
    end

    basefit!(model, constrain=constrain)
    cm = calculation(model)


    if model.variant.adjustment == AC_DXT
        bms_adjust!(model; constrain=constrain && cm == CC_JULIA)
        # bms_adjust!(model; constrain=constrain)
    elseif model.variant.adjustment == AC_E0
        lm_adjust!(model; constrain=constrain && cm == CC_JULIA)
        # lm_adjust!(model; constrain=constrain)
    elseif model.variant.adjustment == AC_DT
        lc_adjust!(model; constrain=constrain && cm == CC_JULIA)
        # lc_adjust!(model; constrain=constrain)
    end

    return model

end

function variation_explained(model::MortalityModel)
    lmxt = logrates(model)
    α = mapslices(lmt -> mean(lmt[.!isinf.(lmt)]), lmxt, dims=2)
    centered = lmxt .- α

    σ = svd(centered).S
    σ² = σ .^ 2

    ∑σ² = sum(σ²)

    pve = cumsum(σ²) ./ ∑σ²

    return pve
end