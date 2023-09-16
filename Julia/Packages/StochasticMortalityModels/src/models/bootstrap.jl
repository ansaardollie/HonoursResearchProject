
include("./../requirements.jl")
include("./../hmd/hmd.jl")
include("./../utilities/utilities.jl")
include("./models.jl")

using ThreadsX
# function deviance_residual(obs, fit)
#     if obs == 0 || fit == 0
#         return 0
#     else
#         return sign(obs - fit) * sqrt(2 * ((obs * log(obs / fit)) - (obs - fit)))
#     end
# end


# function deviance(obs, fit)

#     if obs == 0 || fit == 0
#         return 0
#     else
#         return 2 * (obs * log(obs / fit) - (obs - fit))
#     end
# end


function bootstrap_logrates(model::MortalityModel; n::Int=5000)
    obs = logrates(model, DS_TRAIN)
    fit = mxt_hat(model.parameters; log_scale=true).values

    X = size(obs, 1)
    T = size(obs, 2)

    residuals = obs .- fit

    results = Array{Float64,3}(fill(0.0, X, T, n))
    for i in 1:n
        resample = mapslices(row -> sample(row, T, replace=true), residuals, dims=2)

        results[:, :, i] = fit + resample
    end

    return results
end


BootstrapSampleInfo = NamedTuple{(:lmxt, :constrain),Tuple{Matrix{Float64},Bool}}

function gen_bootstrap_fit(model::MortalityModel)

    Ext = exposures(model, DS_TRAIN)
    dxt = deaths(model, DS_TRAIN)
    ext = deaths(model, DS_TRAIN)
    x = ages(model, DS_TRAIN)
    t = years(model, DS_TRAIN)
    gender = sex(model)
    adjust = adjustment(model)
    cmode = calculation(model)

    if adjust == AC_DT
        gen_func_lc(info::BootstrapSampleInfo) = begin
            lmxt = info.lmxt
            (α′, β′, κ′) = basefit!(lmxt; constrain=info.constrain, cmode=cmode)
            (α, β, κ) = lc_adjust(α′, β′, κ′, t, dxt, Ext; constrain=info.constrain, cmode=cmode)

            return (alphas=α, betas=β, kappas=κ)
        end

        return gen_func_lc
    elseif adjust == AC_DXT
        gen_func_bms(info::BootstrapSampleInfo) = begin
            lmxt = info.lmxt
            (α′, β′, κ′) = basefit!(lmxt; constrain=info.constrain, cmode=cmode)
            (α, β, κ) = bms_adjust(α′, β′, κ′, t, dxt, Ext; constrain=info.constrain, cmode=cmode)
            return (alphas=α, betas=β, kappas=κ)
        end

        return gen_func_bms
    elseif adjust == AC_E0
        gen_func_lm(info::BootstrapSampleInfo) = begin
            lmxt = info.lmxt
            (α′, β′, κ′) = basefit!(lmxt; constrain=info.constrain, cmode=cmode)
            mxt = exp.(lmxt)
            (α, β, κ) = lm_adjust(α′, β′, κ′, t, x, mxt, ext, gender; constrain=info.constrain, cmode=cmode)
            return (alphas=α, betas=β, kappas=κ)
        end

        return gen_func_lm
    else
        gen_func(info::BootstrapSampleInfo) = begin
            lmxt = info.lmxt

            (α, β, κ) = basefit!(lmxt; constrain=info.constrain, cmode=cmode)
            return (alphas=α, betas=β, kappas=κ)
        end

        return gen_func_none
    end

end


function bootstrap(model::MortalityModel; n::Int=5000, constrain::Bool=true)

    lmxt_bsamples = bootstrap_logrates(model, n=n)
    nx = length(ages(model))
    nt = length(years(model))
    data = Dict{Int,BootstrapSampleInfo}()

    for z in axes(lmxt_bsamples, 3)
        lmxt_z = lmxt_bsamples[:, :, z]

        info::BootstrapSampleInfo = (lmxt=lmxt_z, constrain=constrain)

        data[z] = info
    end

    est_func = gen_bootstrap_fit(model)

    output = ThreadsX.map(data) do pair
        return pair[1] => est_func(pair[2])
    end

    sort!(output, by=p -> p[1])

    alpha_samples = map(x -> x[2].alphas, output)
    beta_samples = map(x -> x[2].betas, output)
    kappa_samples = map(x -> x[2].kappas, output)

    α = reshape(collect(Iterators.flatten(alpha_samples)), nx, n)
    β = reshape(collect(Iterators.flatten(beta_samples)), nx, n)
    κ = reshape(collect(Iterators.flatten(kappa_samples)), nt, n)

    return (alphas=α, betas=β, kappas=κ)
    μ_α = mean(α, dims=2)
    μ_β = mean(β, dims=2)
    μ_κ = mean(κ, dims=2)

    # σ_α = std(α, dims=2)
    # σ_β = std(β, dims=2)
    # σ_κ = std(κ, dims=2)

    # med_α = mapslices(x -> quantile(x, 0.5), α, dims=2)
    # med_β = mapslices(x -> quantile(x, 0.5), β, dims=2)
    # med_κ = mapslices(x -> quantile(x, 0.5), κ, dims=2)
end

