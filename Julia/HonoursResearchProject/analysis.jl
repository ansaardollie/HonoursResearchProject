using StochasticMortalityModels

using Gadfly

using PrettyTables
using Dates
using DataFrames
using DataFramesMeta
using LaTeXStrings

function get_file_path(name::String)
    cd = pwd()
    subfolder = string(Dates.today())
    dirpath = joinpath(cd, "outputs", subfolder)
    mkpath(dirpath)
    filepath = joinpath(dirpath, "$name")
    touch(filepath)
    return filepath
end

sex_color_palette = ["red", "blue"]
model_color_palette = ["#f65e94", "#1baa52", "#ebd918", "#30a0c5"]

figure_theme = Theme(
    panel_fill="white",
    background_color="white",
    plot_padding=[1cm],
    highlight_width=0mm,
    guide_title_position=:left,
    colorkey_swatch_shape=:circle,
    key_position=:top,
)
Gadfly.push_theme(figure_theme)
#====================================================
Start [Part 1]: Determine which countries to use
====================================================#

countries = [
    "Australia",
    "Austria",
    "Belarus",
    "Belgium",
    "Bulgaria",
    "Canada",
    "Chile",
    "Croatia",
    "Czechia",
    "Denmark",
    "Estonia",
    "Finland",
    "France",
    "Germany",
    "Greece",
    "Hong Kong",
    "Hungary",
    "Iceland",
    "Ireland",
    "Israel",
    "Italy",
    "Japan",
    "Latvia",
    "Lithuania",
    "Luxembourg",
    "Netherlands",
    "New Zealand",
    "Norway",
    "Poland",
    "Portugal",
    "South Korea",
    "Russia",
    "Slovakia",
    "Slovenia",
    "Spain",
    "Sweden",
    "Switzerland",
    "Taiwan",
    "Ukraine",
    "UK",
    "USA"]


codes = Dict{String,String}(
    "Australia" => "AUS",
    "Austria" => "AUT",
    "Belarus" => "BLR",
    "Belgium" => "BEL",
    "Bulgaria" => "BGR",
    "Canada" => "CAN",
    "Chile" => "CHL",
    "Croatia" => "HRV",
    "Czechia" => "CZE",
    "Denmark" => "DNK",
    "Estonia" => "EST",
    "Finland" => "FIN",
    "France" => "FRA",
    "Germany" => "DEU",
    "Greece" => "GRC",
    "Hong Kong" => "HKG",
    "Hungary" => "HUN",
    "Iceland" => "ISL",
    "Ireland" => "IRL",
    "Israel" => "ISR",
    "Italy" => "ITA",
    "Japan" => "JPN",
    "Latvia" => "LVA",
    "Lithuania" => "LTU",
    "Luxembourg" => "LUX",
    "Netherlands" => "NLD",
    "New Zealand" => "NZL",
    "Norway" => "NOR",
    "Poland" => "POL",
    "Portugal" => "PRT",
    "Republic of Korea" => "KOR",
    "Russia" => "RUS",
    "Slovakia" => "SVK",
    "Slovenia" => "SVN",
    "Spain" => "ESP",
    "Sweden" => "SWE",
    "Switzerland" => "CHE",
    "Taiwan" => "TWN",
    "Ukraine" => "URK",
    "United Kingdom" => "GBR",
    "United States of America" => "USA"
)

pve_by_country = Dict{String,Vector{Float64}}()
exclusions = Dict{String,Int}()
for c in countries
    println()
    println("============================================")
    println("Analysing ", c)
    b_mod = MortalityModel(c, SEX_BOTH; remove_missing=true)
    m_mod = MortalityModel(c, SEX_MALE; remove_missing=true)
    f_mod = MortalityModel(c, SEX_FEMALE; remove_missing=true)
    nt = length(years(b_mod, DS_COMPLETE))


    if nt < 70
        println("Excluding ", c, " insufficient data")
        exclusions[c] = nt
        println("============================================")
        println()
        continue
    end
    try
        pve_b = variation_explained(b_mod)[1]
        pve_m = variation_explained(m_mod)[1]
        pve_f = variation_explained(f_mod)[1]
        pve_by_country[c] = [pve_b, pve_f, pve_m]
    catch e
        println("Issue for country: ", c)
        println(e)
    end
    println("============================================")
    println()
end

spve = sort([p[1] => pve_by_country[p[1]][1] for p in pve_by_country], by=p -> 1 - p[2])

dm = [pve_by_country[p[1]] for p in spve]
dm = round.(mapreduce(permutedims, vcat, dm) * 100, digits=3)

names = map(p -> p[1], spve)



file = get_file_path("Variation Explained.tex")

open(file, "w") do fio

    myio = IOContext(fio, :file => true)

    conf = latex_table_config(myio;
        headers=[:Both, :Female, :Male],
        title="PVE by first singular value",
        rows=names,
        row_label_title="Country",
        table_format=tf_latex_booktabs,
        row_label_alignment=:l,
        formatters=gen_seperator(3),
        hlines=[:begin, :begin, :header, :end, :end])
    pretty_table(myio, dm; conf...)
end

selected_countries = ["Netherlands", "Spain", "Slovakia", "Hungary"]
start_years_lc = [1925, 1925, 1950, 1950]
start_years_lm = [1950, 1950, 1950, 1950]
start_years_bms = [1950, 1950, 1950, 1950]


model_variants = ["LC" "LM" "BMS₁" "BMS₂"
    ; "LC" "LM" "BMS₁" "BMS₂";;;
    "LC" "LM" "BMS₁" "BMS₂"
    ; "LC" "LM" "BMS₁" "BMS₂";;;
    "LC" "LM" "BMS₁" "BMS₂"
    ; "LC" "LM" "BMS₁" "BMS₂";;;
    "LC" "LM" "BMS₁" "BMS₂";
    "LC" "LM" "BMS₁" "BMS₂"]

model_sex = [
    "Female" "Female" "Female" "Female"
    ; "Male" "Male" "Male" "Male";;;
    "Female" "Female" "Female" "Female"
    ; "Male" "Male" "Male" "Male";;;
    "Female" "Female" "Female" "Female"
    ; "Male" "Male" "Male" "Male";;;
    "Female" "Female" "Female" "Female"
    ; "Male" "Male" "Male" "Male"
]

model_country = [
    "Netherlands" "Netherlands" "Netherlands" "Netherlands"
    ; "Netherlands" "Netherlands" "Netherlands" "Netherlands";;;
    "Spain" "Spain" "Spain" "Spain"
    ; "Spain" "Spain" "Spain" "Spain";;;
    "Slovakia" "Slovakia" "Slovakia" "Slovakia"
    ; "Slovakia" "Slovakia" "Slovakia" "Slovakia";;;
    "Hungary" "Hungary" "Hungary" "Hungary";
    "Hungary" "Hungary" "Hungary" "Hungary"
]
#====================================================
End [Part 1]: Determine which countries to use
====================================================#


#====================================================
Start [Part 2]: Create Models for selected countries
====================================================#

models = Array{MortalityModel,3}(undef, 2, 4, 4)

for sci in eachindex(selected_countries)
    sc = selected_countries[sci]
    sy_lc = start_years_lc[sci]
    sy_lm = start_years_lm[sci]
    sy_bms = start_years_bms[sci]
    m1 = MortalityModel(
        sc,
        SEX_MALE;
        remove_missing=true,
        train_years=sy_lc:1999,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=lc
    )
    m2 = MortalityModel(
        sc,
        SEX_MALE;
        remove_missing=true,
        train_years=sy_lm:1999,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=lm
    )
    m3 = MortalityModel(
        sc,
        SEX_MALE;
        remove_missing=true,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=bms
    )
    choose_period!(m3)
    m4 = MortalityModel(
        sc,
        SEX_MALE;
        remove_missing=true,
        train_years=sy_bms:1999,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=bms
    )
    f1 = MortalityModel(
        sc,
        SEX_FEMALE;
        remove_missing=true,
        train_years=sy_lc:1999,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=lc
    )
    f2 = MortalityModel(
        sc,
        SEX_FEMALE;
        remove_missing=true,
        train_years=sy_lm:1999,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=lm
    )
    f3 = MortalityModel(
        sc,
        SEX_FEMALE;
        remove_missing=true,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=bms
    )
    choose_period!(f3)
    f4 = MortalityModel(
        sc,
        SEX_FEMALE;
        remove_missing=true,
        train_years=sy_bms:1999,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=bms
    )
    country_models = [f1 f2 f3 f4; m1 m2 m3 m4]
    models[:, :, sci] = country_models
end


#====================================================
End [Part 2]: Create Models for selected countries
====================================================#


#====================================================
Start [Part 2]: Fit Models & Create Plots/Tables
====================================================#


for m in models
    try
        fit!(m)


    catch e
        l = location(m)
        s = sex(m)
        v = adjustment(m)

        println("Error while fitting ", l, " ", s, " ", v)
        println(e)
    end
end



age_values = collect(Iterators.flatten([ages(m) for m in models]))
alpha_values = collect(Iterators.flatten([alphas(m) for m in models]))
beta_values = collect(Iterators.flatten([betas(m) for m in models]))

nx = 111

variant_values = collect(Iterators.flatten([repeat([v], nx) for v in model_variants]))
sex_values = collect(Iterators.flatten([repeat([s], nx) for s in model_sex]))
country_values = collect(Iterators.flatten([repeat([c], nx) for c in model_country]))

age_dependent_parameters = DataFrame(
    :Age => age_values,
    :Alpha => alpha_values,
    :Beta => beta_values,
    :Country => country_values,
    :Sex => sex_values,
    :Variant => variant_values
)

x_ticks = collect(0:25:110)
y_ticks = collect(0:-2:-10)
p1v1 = plot(
    age_dependent_parameters,
    x="Age",
    y="Alpha",
    ygroup="Country",
    xgroup="Variant",
    color="Sex",
    Scale.x_continuous(minvalue=0, maxvalue=110),
    Scale.color_discrete_manual(sex_color_palette...),
    Guide.ylabel("Estimated αₓ", orientation=:vertical),
    Guide.xlabel("Age (years)"),
    Guide.colorkey(title="Sex", labels=["Female", "Male"]),
    Geom.subplot_grid(
        Geom.line,
        Guide.xticks(ticks=x_ticks),
        Guide.yticks(ticks=y_ticks)
    ),
    Guide.title("Estimated αₓ By Country, Model & Sex"),
);
fp1 = get_file_path("Alpha Plot [Sex Legend].svg")
g1 = SVG(fp1, 40cm, 40cm)

draw(g1, p1v1);

p1v2 = plot(
    age_dependent_parameters,
    x="Age",
    y="Alpha",
    xgroup="Sex",
    ygroup="Country",
    color="Variant",
    Scale.x_continuous(minvalue=0, maxvalue=110),
    Scale.color_discrete_manual(model_color_palette...),
    Guide.ylabel("Estimated αₓ", orientation=:vertical),
    Guide.xlabel("Age (years)"), Guide.title("Estimated αₓ By Country, Sex & Model"),
    Guide.colorkey(; title="Model"),
    Geom.subplot_grid(
        Geom.line,
        Guide.xticks(ticks=x_ticks),
        Guide.yticks(ticks=y_ticks)
    )
);

fp2 = get_file_path("Alpha Plot [Variant Legend].svg")
g2 = SVG(fp2, 20cm, 40cm)
draw(g2, p1v2)


x_ticks = collect(0:25:110)
y_ticks = sort(cat(collect(-0.08:0.05:0.17), 0, dims=1))
p2v1 = plot(
    age_dependent_parameters,
    x="Age",
    y="Beta",
    ygroup="Country",
    xgroup="Variant",
    color="Sex",
    Scale.x_continuous(minvalue=0, maxvalue=110),
    Scale.color_discrete_manual(sex_color_palette...),
    Guide.ylabel("Estimated βₓ", orientation=:vertical),
    Guide.xlabel("Age (years)"),
    Guide.title("Estimated βₓ By Country, Model Variant & Sex"),
    Geom.subplot_grid(
        Geom.line,
        Guide.xticks(ticks=x_ticks),
        Guide.yticks(ticks=y_ticks)
    )
);


fp3 = get_file_path("Beta Plot [Sex Legend].svg")
g3 = SVG(fp3, 40cm, 40cm)
draw(g3, p2v1)

time_dependent_parameters = DataFrame(
    :Year => Int[],
    :Kappa => Float64[],
    :KappaRaw => Float64[],
    :Country => String[],
    :Sex => String[],
    :Variant => String[]
)

for ci in eachindex(models)
    m = models[ci]
    s = ["Female", "Male"][Int(sex(m))-1]
    c = location(m)
    v = model_variants[ci]
    ts = years(m)
    kt = kappas(m)
    krt = kappas(m, PV_UNADJUSTED)

    for i in eachindex(ts)
        t = ts[i]
        κ = kt[i]
        kr = krt[i]
        push!(time_dependent_parameters, [t, κ, kr, c, s, v])
    end
end




x_ticks = collect(1925:15:2000)
y_ticks = collect(-150:50:100)
p3 = plot(
    time_dependent_parameters,
    x="Year",
    y="Kappa",
    xgroup="Sex",
    ygroup="Country",
    color="Variant",
    Scale.x_continuous(minvalue=1925, maxvalue=2000),
    Scale.color_discrete_manual(model_color_palette...),
    Guide.ylabel("Estimated κₜ", orientation=:vertical),
    Guide.title("Estimated Post-Adjustment κₜ By Sex, Country & Model Variant"),
    Guide.colorkey(; title="Model"),
    # Guide.colorkey(; title="Variant", labels=["A", "B", "C"]),
    Guide.xlabel("Time (year)"),
    Geom.subplot_grid(
        Geom.line,
        Guide.xticks(ticks=x_ticks),
        Guide.yticks(ticks=y_ticks)
    )
);


fp5 = get_file_path("Kappa Plot [Post Adjustment].svg")
g5 = SVG(fp5, 20cm, 40cm)
draw(g5, p3)


x_ticks = collect(1925:15:2000)
y_ticks = collect(-100:50:150)

p4 = plot(
    time_dependent_parameters,
    x="Year",
    y="KappaRaw",
    xgroup="Sex",
    ygroup="Country",
    color="Variant",
    Scale.x_continuous(minvalue=1925, maxvalue=2000),
    Scale.color_discrete_manual(model_color_palette...),
    Guide.ylabel("Estimated κₜ", orientation=:vertical),
    Guide.title("Estimated Pre-Adjustment κₜ By Sex, Country & Model Variant"),
    Guide.colorkey(; title="Model"),
    # Guide.colorkey(; title="Variant", labels=["A", "B", "C"]),
    Guide.xlabel("Time (year)"),
    Geom.subplot_grid(
        Geom.line,
        Guide.xticks(ticks=x_ticks),
        Guide.yticks(ticks=y_ticks)
    )
);


fp6 = get_file_path("Kappa Plot [Pre-Adjustment].svg")
g6 = SVG(fp6, 20cm, 40cm)
draw(g6, p4)


rs_data = DataFrame(
    :Year => Int[],
    :Country => String[],
    :Sex => String[],
    :Rs => Float64[]
)


for sci in eachindex(selected_countries)
    sc = selected_countries[sci]

    m3 = MortalityModel(
        sc,
        SEX_MALE;
        remove_missing=true,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=bms
    )
    rm = choose_period!(m3)

    f3 = MortalityModel(
        sc,
        SEX_FEMALE;
        remove_missing=true,
        test_years=2000:2019,
        calculation_mode=CC_JULIA,
        variant=bms
    )
    rf = choose_period!(f3)
    for i in axes(rm, 1)
        y = rm[i, 1]
        rs = rm[i, 2]
        push!(rs_data, [y, sc, "Male", rs])
    end

    for i in axes(rf, 1)
        y = rf[i, 1]
        rs = rf[i, 2]
        push!(rs_data, [y, sc, "Female", rs])
    end


end


plts = Array{Plot,2}(undef, 4, 2)

sexes = ["Female", "Male"]
for si in eachindex(sexes)
    s = sexes[si]
    for sci in eachindex(selected_countries)
        sc = selected_countries[sci]
        sdf = @subset(rs_data, :Country .== sc, :Sex .== s)
        min_t = sdf.Year[findmin(sdf.Rs)[2]]
        min_rs = findmin(sdf.Rs)[1]
        start_year = Int(floor(minimum(sdf.Year) / 20) * 20)
        end_year = Int(ceil(maximum(sdf.Year) / 20) * 20)

        main_ticks = collect(start_year:20:end_year)
        issue_idx = findall(x -> abs(x - min_t) < 5, main_ticks)
        deleteat!(main_ticks, issue_idx)
        xt = sort!([min_t, 2000, main_ticks...])

        println("[$sc / $s] years = $(xt[1]) -> $(xt[end])")
        l1 = layer(
            sdf,
            x=:Year,
            y=:Rs,
            color=[s == "Male" ? colorant"blue" : colorant"red"],
            Geom.line,
        )

        l2 = layer(
            x=[min_t],
            y=[min_rs],
            color=[colorant"black"],
            shape=["cross"]
        )

        l3 = layer(
            x=[min_t, start_year],
            y=[0, min_rs],
            xend=[min_t, min_t],
            yend=[min_rs, min_rs],
            color=[colorant"green"],
            Geom.segment
        )

        plt = plot(l1, l2, l3,
            Guide.xlabel("Starting Year (S)", orientation=:horizontal),
            Guide.ylabel("R(S)"),
            Guide.xticks(ticks=xt, orientation=:horizontal),
            Scale.x_continuous(minvalue=start_year, maxvalue=2000),
            style(
                key_position=:none,
                minor_label_font_size=5pt
            ),
            Guide.title("$sc / $s"))

        plts[sci, si] = plt
    end
end

bms_plot = title(gridstack(plts), "Mean Deviance Ratio Statistics")

fp_bms = get_file_path("BMS Selection Grouped.svg")
g_bms = SVG(fp_bms, 20cm, 40cm)

draw(g_bms, bms_plot)

#====================================================
End [Part 2]: Fit Models & Create Plots/Tables
====================================================#


#====================================================
Start [Part 3]: Forecast Models & Plot
====================================================#

forecasts_io = Array{ModelForecasts,3}(undef, 2, 4, 4)

for mi in eachindex(models)
    m = models[mi]
    fm = forecast(m; confidence_level=0.95, uncertainty=UM_INNOVATION_ONLY)
    forecasts_io[mi] = fm
end

forecasts_id = Array{ModelForecasts,3}(undef, 2, 4, 4)

for mi in eachindex(models)
    m = models[mi]
    fm = forecast(m; confidence_level=0.95, uncertainty=UM_INNOVDRIFT)
    forecasts_id[mi] = fm
end


le_forecast_data2 = DataFrame(
    :Year => Int[],
    :Country => String[],
    :Sex => String[],
    :Variant => String[],
    :TrueValue => Float64[],
    :Forecast => Float64[],
    :UpperBound1 => Float64[],
    :LowerBound1 => Float64[],
    :BoundType => String[],
)

for idx in eachindex(models)
    m = models[idx]
    fd_io = forecasts_io[idx]
    fd_id = forecasts_id[idx]


    c = location(m)
    s = ["Female", "Male"][Int(sex(m))-1]
    v = model_variants[idx]

    train_years = years(m, DS_TRAIN)[end-4:end]
    test_years = years(m, DS_TEST)


    train_e0 = lifespans(m, DS_TRAIN)[1, end-4:end]
    test_e0 = lifespans(m, DS_TEST)[1, :]
    forecasted_e0 = fd_io.lifespans.forecast.values[1, :]
    ub1_e0 = fd_io.lifespans.upper.values[1, :]
    lb1_e0 = fd_io.lifespans.lower.values[1, :]
    forecasted_e0 = fd_id.lifespans.forecast.values[1, :]
    ub2_e0 = fd_id.lifespans.upper.values[1, :]
    lb2_e0 = fd_id.lifespans.lower.values[1, :]

    for i in eachindex(test_years)
        t = test_years[i]
        true_e0 = test_e0[i]
        fc = forecasted_e0[i]
        ub1 = ub1_e0[i]
        lb1 = lb1_e0[i]

        push!(le_forecast_data2, [t, c, s, v, true_e0, fc, ub1, lb1, "Innovation Only"])
    end


    for i in eachindex(test_years)
        t = test_years[i]
        true_e0 = test_e0[i]
        fc = forecasted_e0[i]
        ub1 = ub2_e0[i]
        lb1 = lb2_e0[i]

        push!(le_forecast_data2, [t, c, s, v, true_e0, fc, ub1, lb1, "Innovation + Drift"])
    end


end



for sc in selected_countries
    myticks = sc == "Netherlands" ?
              collect(66:2:82) :
              sc == "Spain" ?
              collect(74:2:84) :
              nothing
    fyticks = sc == "Netherlands" ?
              collect(72:2:88) :
              sc == "Spain" ?
              collect(82:1:88) :
              sc == "Slovakia" ?
              collect(76:1:84) :
              collect(74:1:82)
    plt_stacks = Matrix{Plot}(undef, 2, 4)
    for s in ["Male", "Female"]
        sidx = s == "Male" ? 2 : 1
        yticks = s == "Male" ? myticks : fyticks
        for v in ["LC", "LM", "BMS₁", "BMS₂"]
            vidx = v == "LC" ?
                   1 :
                   v == "LM" ?
                   2 :
                   v == "BMS₁" ?
                   3 : 4
            varname = vidx == 3 ? "BMS₁" :
                      vidx == 4 ? "BMS₂" :
                      v
            subdf = @subset(le_forecast_data2, :Country .== sc, :Sex .== s, :Variant .== v)

            fc1_layer = layer(
                x="Year",
                y="Forecast",
                ymin="LowerBound1",
                ymax="UpperBound1",
                color="BoundType",
                alpha="BoundType",
                Geom.line,
                Geom.ribbon(; fill=true)
            )

            tv_layer = layer(
                x="Year",
                y="TrueValue",
                color=[colorant"red"],
                size=[0.5mm],
                Geom.point()
            )

            plt = !isnothing(yticks) ?
                  plot(
                subdf,
                tv_layer,
                fc1_layer,
                Guide.xlabel("Time (year)"),
                Guide.ylabel("e₀ Forecast"),
                Guide.title("$varname for $s"),
                Scale.color_discrete_manual("blue", "green2"),
                Guide.yticks(ticks=yticks),
                style(
                    alphas=[1, 0.5],
                    key_position=:none
                )
            ) :
                  plot(
                subdf,
                tv_layer,
                fc1_layer,
                Guide.xlabel("Time (year)"),
                Guide.ylabel("e₀ Forecast"),
                Guide.title("$varname for $s"),
                Scale.color_discrete_manual("blue", "green2"),
                style(
                    alphas=[1, 0.5],
                    key_position=:none
                )
            )
            plt_stacks[sidx, vidx] = plt

        end

    end

    country_plt = title(gridstack(plt_stacks), "Forecasted Life Expectancy At Birth for $sc")

    f = get_file_path("Life Expectancy Forecasts [$sc].svg")
    g = SVG(f, 40cm, 20cm)
    draw(g, country_plt)

end


#====================================================
End [Part 3]: Forecast Models & Plot
====================================================#


#====================================================
Start [Part 4]: Forecast Accuracy
====================================================#

AccuracyMeasure = Dict{String,Vector{Float64}}


accuracies = Array{AccuracyMeasure,3}(undef, 2, 4, 4)

for mi in eachindex(models)
    m = models[mi]
    fd = forecasts_io[mi]

    mase_r = maswe(m, fd; max_lookahead=25, variable=FV_RATE)
    mase_mrl = maswe(m, fd; max_lookahead=25, variable=FV_MRL)
    mase_leab = maswe(m, fd; max_lookahead=25, variable=FV_LEAB)
    am::AccuracyMeasure = Dict{String,Vector{Float64}}(
        "Rates" => mase_r,
        "MeanResidualLifetimes" => mase_mrl,
        "LEAB" => mase_leab
    )

    accuracies[mi] = am
end

accuracy_data = DataFrame(
    :Horizon => Int[],
    :Country => String[],
    :Sex => String[],
    :Variant => String[],
    :RMase => Float64[],
    :MrlMase => Float64[],
    :LeabMase => Float64[]
)


for mi in eachindex(models)
    m = models[mi]
    c = location(m)
    s = ["Female", "Male"][Int(sex(m))-1]
    c = location(m)
    v = model_variants[mi]

    ram = accuracies[mi]["Rates"]
    mrlam = accuracies[mi]["MeanResidualLifetimes"]
    leabam = accuracies[mi]["LEAB"]
    for i in eachindex(lram)
        ry = ram[i]
        mrly = mrlam[i]
        leaby = leabam[i]
        push!(accuracy_data, [i, c, s, v, ry, mrly, leaby])
    end
end

r_fa_plt = plot(
    accuracy_data,
    x="Horizon",
    y="RMase",
    color="Variant",
    ygroup="Country",
    xgroup="Sex",
    yintercept=[1],
    size=[0.5mm],
    Guide.title("Forecast Accuracy of Mortality Rates"),
    Scale.x_continuous(minvalue=0, maxvalue=25),
    Scale.color_discrete_manual(model_color_palette...),
    Guide.ylabel("MASWEₖ", orientation=:vertical),
    Guide.xlabel("Naive Forecast Step (k)"),
    Geom.subplot_grid(
        Geom.hline,
        Geom.point,
        Guide.xticks(ticks=x_ticks),
        free_y_axis=true
    )
);

fp15 = get_file_path("Accuracy Measures [Rates].svg")
g15 = SVG(fp15, 20cm, 40cm)
draw(g15, r_fa_plt)

le_fa_plt = plot(
    accuracy_data,
    x="Horizon",
    y="MrlMase",
    color="Variant",
    ygroup="Country",
    xgroup="Sex",
    yintercept=[1],
    size=[0.5mm],
    Scale.x_continuous(minvalue=0, maxvalue=25),
    Scale.color_discrete_manual(model_color_palette...),
    Guide.title("Forecast Accuracy of Expected Residual Lifetime"),
    Guide.ylabel("MASWEₖ", orientation=:vertical),
    Guide.xlabel("Naive Forecast Step (k)"),
    Geom.subplot_grid(
        Geom.hline,
        Geom.point,
        Guide.xticks(ticks=x_ticks),
        free_y_axis=true
    )
);

fp16 = get_file_path("Accuracy Measures [Mean Residual Lifetimes].svg")
g16 = SVG(fp16, 20cm, 40cm)
draw(g16, le_fa_plt)

leab_fa_plt = plot(
    accuracy_data,
    x="Horizon",
    y="LeabMase",
    color="Variant",
    ygroup="Country",
    xgroup="Sex",
    yintercept=[1],
    size=[0.5mm],
    Guide.title("Forecast Accuracy of Life Expectancy At Birth"),
    Scale.x_continuous(minvalue=0, maxvalue=25),
    Scale.color_discrete_manual(model_color_palette...),
    Guide.ylabel("MASEₖ", orientation=:vertical),
    Guide.xlabel("Naive Forecast Step (k)"),
    Geom.subplot_grid(
        Geom.hline,
        Geom.point,
        Guide.xticks(ticks=x_ticks),
        free_y_axis=true
    )
);


fp18 = get_file_path("Accuracy Measures [LEAB].svg")
g18 = SVG(fp18, 20cm, 40cm)
draw(g18, leab_fa_plt)


r_info = Matrix{Float64}(undef, 8, 8)
mrl_info = Matrix{Float64}(undef, 8, 8)
leab_info = Matrix{Float64}(undef, 8, 8)

idx = [(1, 1) (1, 2) (1, 3) (1, 4)
    ; (5, 1) (5, 2) (5, 3) (5, 4);;;
    (2, 1) (2, 2) (2, 3) (2, 4)
    ; (6, 1) (6, 2) (6, 3) (6, 4);;;
    (3, 1) (3, 2) (3, 3) (3, 4)
    ; (7, 1) (7, 2) (7, 3) (7, 4);;;
    (4, 1) (4, 2) (4, 3) (4, 4);
    (8, 1) (8, 2) (8, 3) (8, 4)]


for mi in eachindex(models)
    m = models[mi]
    c = location(m)
    s = ["Female", "Male"][Int(sex(m))-1]
    c = location(m)
    v = model_variants[mi]

    ram = accuracies[mi]["Rates"]
    ram_bp = findfirst(a -> a < 1, ram)
    ram_bp = isnothing(ram_bp) ? 0 : ram_bp

    leam = accuracies[mi]["MeanResidualLifetimes"]
    leam_bp = findfirst(a -> a < 1, leam)
    leam_bp = isnothing(leam_bp) ? 0 : leam_bp

    leab = accuracies[mi]["LEAB"]
    leab_bp = findfirst(a -> a < 1, leab)
    leab_bp = isnothing(leab_bp) ? 0 : leab_bp

    ri = idx[mi][1]
    ci = idx[mi][2]
    cibp = ci + 4

    r_info[ri, ci] = ram[1]
    r_info[ri, cibp] = ram_bp

    mrl_info[ri, ci] = leam[1]
    mrl_info[ri, cibp] = leam_bp

    leab_info[ri, ci] = leab[1]
    leab_info[ri, cibp] = leab_bp

end


am_info_names = ["Netherlands", "Spain", "Slovakia", "Hungary", "Netherlands", "Spain", "Slovakia", "Hungary"]

file = get_file_path("Rates MASE.tex")

sep = gen_seperator(3)

ri_dm = [lri == 0 ? L">25" : isinteger(lri) ? L"%$(Int(lri))" : L"%$(sep(lri, 1, 1))" for lri in r_info]
open(file, "w") do fio

    myio = IOContext(fio, :file => true)

    conf = latex_table_config(myio;
        headers=[:LC, :LM, :BMS1, :BMS2, :LC, :LM, :BMS1, :BMS2],
        title="Rates Forecast Accuracy per Model,Country,Sex",
        rows=am_info_names,
        row_label_title="Country",
        table_format=tf_latex_booktabs,
        row_label_alignment=:l,
        type=:longtable,
        hlines=[:begin, :begin, :header, :end, :end])
    pretty_table(myio, ri_dm; conf...)
end


file = get_file_path("MRL MASE.tex")

sep = gen_seperator(3)

mrl_dm = [lri == 0 ? L">25" : isinteger(lri) ? L"%$(Int(lri))" : L"%$(sep(lri, 1, 1))" for lri in mrl_info]
open(file, "w") do fio

    myio = IOContext(fio, :file => true)

    conf = latex_table_config(myio;
        headers=[:LC, :LM, :BMS1, :BMS2, :LC, :LM, :BMS1, :BMS2],
        title="Expected Residual Lifetime Forecast Accuracy per Model,Country,Sex",
        rows=am_info_names,
        row_label_title="Country",
        table_format=tf_latex_booktabs,
        row_label_alignment=:l,
        type=:longtable,
        hlines=[:begin, :begin, :header, :end, :end])
    pretty_table(myio, mrl_dm; conf...)
end


file = get_file_path("LEAB MASE.tex")

sep = gen_seperator(3)

leab_dm = [lri == 0 ? L">25" : isinteger(lri) ? L"%$(Int(lri))" : L"%$(sep(lri, 1, 1))" for lri in leab_info]
open(file, "w") do fio

    myio = IOContext(fio, :file => true)

    conf = latex_table_config(myio;
        headers=[:LC, :LM, :BMS1, :BMS2, :LC, :LM, :BMS1, :BMS2],
        title="Life Expectancy At Birth Forecast Accuracy per Model,Country,Sex",
        rows=am_info_names,
        row_label_title="Country",
        table_format=tf_latex_booktabs,
        row_label_alignment=:l,
        type=:longtable,
        hlines=[:begin, :begin, :header, :end, :end])
    pretty_table(myio, leab_dm; conf...)
end


residual_data = DataFrame(
    :Year => Int[],
    :Age => Int[],
    :Country => String[],
    :Sex => String[],
    :Variant => String[],
    :Residual => Float64[],
    :ResidualType => String[]
)

for mi in eachindex(models)
    m = models[mi]

    c = location(m)
    s = ["Female", "Male"][Int(sex(m))-1]
    c = location(m)
    v = model_variants[mi]

    x_ages = ages(m)
    t_years = years(m)

    res_lr = residuals(m, variable=FV_LOGRATE)
    res_r = residuals(m, variable=FV_RATE)
    res_mrl = residuals(m, variable=FV_MRL)
    res_deaths = residuals(m, variable=FV_DEATHS)
    for x in x_ages
        xi = indexin(x, x_ages)[1]
        for y in t_years
            yi = indexin(y, t_years)[1]
            lr = res_lr[xi, yi]
            r = res_r[xi, yi]
            mrl = res_mrl[xi, yi]
            d = res_deaths[xi, yi]
            # @show (c, s, v, xi, yi, lr, r, mrl, d)
            push!(residual_data, [y, x, c, s, v, lr, "LogRate"])
            push!(residual_data, [y, x, c, s, v, r, "Rate"])
            push!(residual_data, [y, x, c, s, v, mrl, "MRL"])
            push!(residual_data, [y, x, c, s, v, d, "Deaths"])
        end
    end

end

subdf = @subset(residual_data, :Country .== "Netherlands", :Sex .== "Male", :Variant .== "LC", :ResidualType .== "Rate")
yticks = collect(0:10:110)
xticks = collect(1925:5:2000)

variants = ["LC", "LM", "BMS₁", "BMS₂"]
sexes = ["Male", "Female"]
variables = ["LogRate", "Rate", "MRL", "Deaths"]
varnames = ["Log Mortality Rate", "Mortality Rate", "Expected Residual Lifetime", "Death Count"]

for si in eachindex(sexes)
    s = sexes[si]
    for vari in eachindex(variables)
        var = variables[vari]
        varn = varnames[vari]

        plots = Matrix{Plot}(undef, 4, 4)
        for sci in eachindex(selected_countries)
            sc = selected_countries[sci]
            for vi in eachindex(variants)
                v = variants[vi]
                subdf = @subset(residual_data, :Country .== sc, :Sex .== s, :Variant .== v, :ResidualType .== var)

                min_y = minimum(subdf.Year)
                max_y = maximum(subdf.Year)
                if (max_y - min_y) < 20
                    rd = 5
                else
                    rd = 10
                end
                start_year = Int(ceil(min_y / rd) * rd)
                end_year = 2000

                xticks = [minimum(subdf.Year), collect(start_year:rd:end_year)...]
                yticks = collect(0:10:110)
                modv = v == "BMS1" ? "BMS₁" : v == "BMS2" ? "BMS₂" : v
                plt = plot(
                    subdf,
                    x=:Year,
                    y=:Age,
                    color=:Residual,
                    Geom.rectbin,
                    Guide.xticks(ticks=xticks, orientation=:vertical),
                    Guide.yticks(ticks=yticks),
                    Guide.xlabel("Time (year)"),
                    Guide.ylabel("Age (years)"),
                    Guide.title("$sc / $modv"),
                    style(
                        key_position=:right
                    )
                )

                plots[sci, vi] = plt

            end
        end

        p1 = title(gridstack(plots), "$varn Residuals for $(s)s")
        fp_plt = get_file_path("Residual Heatmaps [$s - $varn].svg")
        g = SVG(fp_plt, 40cm, 40cm)
        draw(g, p1)
    end
end




forecast_errors_data = DataFrame(
    :Year => Int[],
    :Age => Int[],
    :Country => String[],
    :Sex => String[],
    :Variant => String[],
    :ForecastError => Float64[],
    :Type => String[]
)

for mi in eachindex(models)
    m = models[mi]
    fd = forecasts_io[mi]

    c = location(m)
    s = ["Female", "Male"][Int(sex(m))-1]
    c = location(m)
    v = model_variants[mi]

    x_ages = ages(m, DS_TEST)
    t_years = years(m, DS_TEST)

    res_lr = forecast_errors(m, fd, variable=FV_LOGRATE)
    res_r = forecast_errors(m, fd, variable=FV_RATE)
    res_mrl = forecast_errors(m, fd, variable=FV_MRL)
    res_deaths = forecast_errors(m, fd, variable=FV_DEATHS)
    for x in x_ages
        xi = indexin(x, x_ages)[1]
        for y in t_years
            yi = indexin(y, t_years)[1]
            lr = res_lr[xi, yi]
            r = res_r[xi, yi]
            mrl = res_mrl[xi, yi]
            d = res_deaths[xi, yi]
            # @show (c, s, v, xi, yi, lr, r, mrl, d)
            push!(forecast_errors_data, [y, x, c, s, v, lr, "LogRate"])
            push!(forecast_errors_data, [y, x, c, s, v, r, "Rate"])
            push!(forecast_errors_data, [y, x, c, s, v, mrl, "MRL"])
            push!(forecast_errors_data, [y, x, c, s, v, d, "Deaths"])
        end
    end

end


variants = ["LC", "LM", "BMS₁", "BMS₂"]
sexes = ["Male", "Female"]
variables = ["LogRate", "Rate", "MRL", "Deaths"]
varnames = ["Log Mortality Rate", "Mortality Rate", "Expected Residual Lifetime", "Death Count"]

f_mrl_min = -3
f_mrl_max = 3
m_mrl_min = -10
m_mrl_max = 10
mxt_min = -0.2
mxt_max = 0.2

for si in eachindex(sexes)
    s = sexes[si]
    mrl_min = s == "Male" ? m_mrl_min : f_mrl_min
    mrl_max = s == "Male" ? m_mrl_max : f_mrl_max
    for vari in eachindex(variables)
        var = variables[vari]
        varn = varnames[vari]
        low = vari == 2 ? mxt_min : vari == 3 ? mrl_min : nothing
        high = vari == 2 ? mxt_max : vari == 3 ? mrl_max : nothing
        plots = Matrix{Plot}(undef, 4, 4)
        for sci in eachindex(selected_countries)
            sc = selected_countries[sci]
            for vi in eachindex(variants)
                v = variants[vi]
                subdf = @subset(forecast_errors_data, :Country .== sc, :Sex .== s, :Variant .== v, :Type .== var)

                xticks = collect(2000:2:2020)
                yticks = collect(0:10:110)
                modv = v == "BMS1" ? "BMS₁" : v == "BMS2" ? "BMS₂" : v
                plt = vari in [2, 3] ? plot(
                    subdf,
                    x=:Year,
                    y=:Age,
                    color=:ForecastError,
                    Geom.rectbin,
                    Guide.xticks(ticks=xticks, orientation=:vertical),
                    Guide.yticks(ticks=yticks),
                    Guide.xlabel("Time (year)"),
                    Guide.ylabel("Age (years)"),
                    Guide.title("$sc / $modv"),
                    Guide.colorkey(; title="Error"),
                    Scale.color_continuous(minvalue=low, maxvalue=high),
                    style(
                        key_position=:right
                    )
                ) : plot(
                    subdf,
                    x=:Year,
                    y=:Age,
                    color=:ForecastError,
                    Geom.rectbin,
                    Guide.xticks(ticks=xticks, orientation=:vertical),
                    Guide.yticks(ticks=yticks),
                    Guide.xlabel("Time (year)"),
                    Guide.ylabel("Age (years)"),
                    Guide.title("$sc / $modv"),
                    Guide.colorkey(; title="Error"),
                    style(
                        key_position=:right
                    )
                )


                plots[sci, vi] = plt

            end
        end
        p1 = title(gridstack(plots), "$varn Forecast Errors for $(s)s")
        fp_plt = get_file_path("FE Heatmaps [$s - $varn].svg")
        g = SVG(fp_plt, 40cm, 40cm)
        draw(g, p1)
    end
end


t_years = 1925:1999
x_ages = 0:110

nt = length(t_years)
nx = length(x_ages)

random = randn(nx, nt)

xvals = Vector{Int}()
yvals = Vector{Int}()
resvals = Vector{Float64}()

for xi in eachindex(t_years)
    for yi in eachindex(x_ages)
        age = x_ages[yi]
        year = t_years[xi]
        res = random[yi, xi]

        push!(xvals, year)
        push!(yvals, age)
        push!(resvals, res)

    end
end

rand_plt = plot(
    x=xvals,
    y=yvals,
    color=resvals,
    Geom.rectbin,
    Guide.xticks(ticks=collect(1925:15:2000), orientation=:vertical),
    Guide.yticks(ticks=collect(0:10:110)),
    Guide.xlabel("Time (year)"),
    Guide.ylabel("Age (years)"),
    Guide.title("Random Noise"),
    style(
        key_position=:right
    )
)

# draw(SVG(20cm, 20cm), rand_plt)

fp_rand = get_file_path("Random Noise.svg")
g_rand = SVG(fp_rand, 10cm, 10cm)
draw(g_rand, rand_plt)