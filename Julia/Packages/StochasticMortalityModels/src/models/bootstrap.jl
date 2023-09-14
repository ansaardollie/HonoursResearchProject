obs = rm.deaths.fit.data
ed = fitted_deaths(rm.parameters, rm.exposures.fit).data

include("./models.jl")
using NLsolve
using Roots
using SparseArrays
using Distributions

function deviance_residual(obs, fit)
    if obs == 0 || fit == 0
        return 0
    else
        return sign(obs - fit) * sqrt(2 * ((obs * log(obs / fit)) - (obs - fit)))
    end
end


function deviance(obs, fit)

    if obs == 0 || fit == 0
        return 0
    else
        return 2 * (obs * log(obs / fit) - (obs - fit))
    end
end


function bootstrap_logrates(model::MortalityModel; n::Int=5000)
    obs = model.logrates.fit.data
    fit = mxt_hat(model.parameters; log_scale=true).data

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

# function bootstrap_deaths(obs, fit; n::Int=5000)

#     directions = sign.(obs - fit)
#     residuals = deviance_residual.(obs, fit)
#     X = size(obs, 1)
#     T = size(obs, 2)


#     idx_m = repeat(reshape(1:T, 1, T), X)

#     results = Array{Float64,3}(fill(0, X, T, n))
#     for i in 1:n
#         res_idx = mapslices(row -> sample(row, T, replace=true), idx_m, dims=2)
#         resampled_residuals = [residuals[ci.I[1], res_idx[ci]] for ci in CartesianIndices(res_idx)]
#         resampled_directions = [directions[ci.I[1], res_idx[ci]] for ci in CartesianIndices(res_idx)]


#         w_mtx = (resampled_residuals .^ 2 ./ (2 .* fit)) .- 1
#         w_mtx = [isnan(w) || isinf(w) ? 0 : w for w in w_mtx]

#         ws = vec(w_mtx)

#         x_start = vec(map(CartesianIndices(resampled_directions)) do ci
#             if w_mtx[ci] == -1
#                 return 1
#             elseif resampled_directions[ci] == 0 || w_mtx[ci] == 0
#                 return 0
#             elseif resampled_directions[ci] < 0
#                 return rand(Uniform(0, -w_mtx[ci]))
#             else
#                 return rand(Uniform(-w_mtx[ci], 10 + (-w_mtx[ci])))
#             end
#         end)

#         x_sols = fill(0.0, X * T)
#         for j in eachindex(ws)
#             try
#                 w = ws[j]
#                 f(x) = begin
#                     if w == 0
#                         return 0
#                     else
#                         return x * log(x) - x - w
#                     end
#                 end

#                 f′(x) = begin
#                     if w == 0
#                         return 0
#                     else
#                         return log(x)
#                     end
#                 end

#                 x0 = x_start[j]
#                 if w == 0
#                     x = 0
#                 elseif w == -1
#                     x = 1
#                 else
#                     x = find_zero((f, f′), x0, Roots.Newton())
#                 end
#                 x_sols[j] = x
#             catch e
#                 println(e)
#                 println("Error occured for")
#                 println((Index=j, Intercept=w, X_Start=x0))

#             end
#         end
#         # solve_func(F, J, x_params) = begin
#         #     lx = [ws[i] == 0 ? 0 : log(x_params[i]) for i in eachindex(x_params)]

#         #     fv = [ws[i] == 0 ? 0 : x_params[i] * lx[i] - x_params[i] - ws[i] for i in eachindex(x_params)]

#         #     jd = [ws[i] == 0 ? 0 : lx[i] for i in eachindex(x_params)]

#         #     if !isnothing(F)
#         #         copyto!(F, fv)
#         #     end

#         #     if !isnothing(J)
#         #         for i in 1:size(J, 1), j in 1:size(J, 2)
#         #             J[i, j] = i == j ? jd[i] : 0
#         #         end
#         #     end
#         # end

#         # f!(F, x_params) = begin
#         #     lx = [ws[i] == 0 ? 0 : log(x_params[i]) for i in eachindex(x_params)]

#         #     fv = [ws[i] == 0 ? 0 : x_params[i] * lx[i] - x_params[i] - ws[i] for i in eachindex(x_params)]

#         #     copyto!(F, fv)
#         #     return F
#         # end

#         # j!(J, x_params) = begin
#         #     lx = [ws[i] == 0 ? 0 : log(x_params[i]) for i in eachindex(x_params)]
#         #     jd = [ws[i] == 0 ? 0 : lx[i] for i in eachindex(x_params)]

#         #     r = size(J, 1)
#         #     for k in 1:r
#         #         J[k, k] = jd[k]
#         #     end
#         #     return J
#         # end

#         # F0 = f!(fill(0.0, X * T), x_start)
#         # J0 = j!(spdiagm(fill(0.0, X * T)), x_start)

#         # df = OnceDifferentiable(f!, j!, x_start, F0, J0)

#         # results = nlsolve(df, x_start, method=:newton)
#         bs_deaths = reshape(results.zero, X, T) .* fit
#         results[:, :, i] = bs_deaths
#     end
#     # ov = vec(obs)
#     # fv = vec(fit)

#     # # issue_idx = sort(collect(Set([findall(x -> x == 0, ov)..., findall(x -> x == 0, fv)...])))
#     # resids = deviance_residual.(obs, fit)
#     # init_x = obs ./ fit

#     # init_x[isnan.(init_x)] .= zero(Float64)
#     # X = size(obs, 1)
#     # T = size(obs, 2)
#     # results = Array{Float64,3}(fill(0, X, T, n))

#     # for i in 1:n
#     #     res_i = mapslices(row -> sample(row, T, replace=true), resids, dims=2)
#     #     w_i = (res_i .^ 2 ./ (2 * fit)) .- 1

#     #     w_v = vec(w_i)
#     #     w_v[isnan.(w_v)] .= 0

#     #     bdx = Vector{Float64}(undef, X * T)

#     #     for i in eachindex(bdx)
#     #         w = w_v[i]
#     #         x0 = init_x[i]

#     #         if w == 0
#     #             bdx[i] = 0
#     #             continue
#     #         end

#     #         f(x) = x * log(x) - x - w
#     #         f′(x) = log(x)

#     #         x_solve = find_zero((f, f′), x0, Roots.Newton())
#     #         bdx[i] = x_solve
#     #     end

#     #     # solve_func(F, J, x) = begin

#     #     #     nan_idx = findall(j -> isnan(j), x)
#     #     #     println("Nan indices = ", nan_idx)

#     #     #     inf_idx = findall(j -> isinf(j), x)
#     #     #     println("Inf indices = ", inf_idx)

#     #     #     lx = [y == 0 ? 0 : log(y) for y in x]

#     #     #     inf_lx_idx = findall(j -> isinf(j), lx)
#     #     #     lx[inf_lx_idx] .= 1

#     #     #     if !isnothing(F)
#     #     #         fv = x .* lx .- x .- w_v
#     #     #         nan_f_idx = findall(j -> isnan(j), fv)
#     #     #         println("Nan indices = ", nan_f_idx)

#     #     #         inf_f_idx = findall(j -> isinf(j), fv)
#     #     #         println("Inf indices = ", inf_f_idx)

#     #     #         fv[nan_f_idx] .= 0
#     #     #         fv[inf_f_idx] .= 0

#     #     #         copyto!(F, fv)
#     #     #     end

#     #     #     if !isnothing(J)
#     #     #         J[diagind(J)] = lx
#     #     #     end
#     #     # end

#     #     # results = nlsolve(only_fj!(solve_func), init_x)
#     #     bdx = bdx .* fv
#     #     results[:, :, i] = reshape(bdx, X, T)

#     # end


# end