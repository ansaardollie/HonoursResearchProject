export expected_lifetime, lexpectancies

function expected_lifetime(
    mx, ages;
    sex::Sex,
    at_age=[0],
    mode::CalculationChoice=CC_JULIA,
    radix=1_000_000,
    debug=false)


    function a0(m0, sex::Sex, mode::CalculationChoice)
        config = a0_config(mode, sex)

        (lb, slope, intercept) = config

        idx = mode == CC_JULIA ? 3 : 2
        for i in eachindex(lb)
            if m0 >= lb[i]
                idx = i
                break
            end
        end

        return intercept[idx] + slope[idx] * m0
    end

    n = length(mx)

    ax = Vector{Any}(fill(0.5, n))

    if ages[1] == 0
        m0 = mx[1]
        a0 = a0(m0, sex, mode)
        ax[1] = a0
    end


    if ages[n] == 110
        aend = 1 / mx[n]
        ax[n] = aend
    end

    qx = mx ./ (1 .+ (1 .- ax) .* mx)

    if ages[n] == 110
        qx[n] = 1
    end

    px = 1 .- qx

    lx = cumprod(cat([1], px[1:end-1], dims=1)) .* radix
    dx = lx .* qx

    Lx = lx .- (1 .- ax) .* dx
    Tx = reverse(cumsum(reverse(Lx)))

    ex = Tx ./ lx

    if debug
        df = DataFrame(:Age => ages, :mx => mx, :ax => ax, :qx => qx, :ex => ex)

        println(df)
    end
    return ex[indexin(at_age, ages)]

end


function lexpectancies(
    mxt::Matrix{Float64},
    ages,
    years;
    sex,
    at_age=[0],
    mode::CalculationChoice=CC_JULIA,
    radix=1_000_000,
    debug=false
)

    output = Matrix{Float64}(undef, length(at_age), length(years))

    for i in eachindex(years)
        year = years[i]
        mx = mxt[:, i]

        les = expected_lifetime(mx, ages; sex=sex, mode=mode, at_age=at_age, radix=radix, debug=debug)
        output[:, i] = les
    end

    return output

end