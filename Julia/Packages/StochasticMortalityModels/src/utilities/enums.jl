export Sex, SEX_BOTH, SEX_FEMALE, SEX_MALE
export Variant, VARIANT_BMS, VARIANT_LC, VARIANT_LM
using Core: Argument
export CalculationMode, CM_DEMOGRAPHY, CM_JULIA
export DataDimension, DD_AGE, DD_YEAR
export MortalityDataSource, MDS_OBSERVED, MDS_FITTED, MDS_PREDICTED, MDS_CALCULATED
export MortalityDataCategory, MDC_EXPOSURES, MDC_DEATHS, MDC_APPROXIMATE_DEATHS, MDC_RATES, MDC_LOGRATES, MDC_LIFE_EXPECTANCIES, MDC_PARAMETERS, MDC_OTHER
export sexmatch, mdc_convert, mdc_label, mdc_shortlabel
export a0_config
export UncertaintyMode, UM_INNOVATION_ONLY, UM_INNOVDRIFT
export AdjustmentChoice, AC_NONE, AC_DT, AC_E0, AC_DXT
export JumpOffRate, JR_FITTED, JF_ACTUAL

@enum Sex begin
    SEX_BOTH = 1
    SEX_FEMALE = 2
    SEX_MALE = 3
end

@enum Variant begin
    VARIANT_BMS = 1
    VARIANT_LC = 2
    VARIANT_LM = 3
end

@enum AdjustmentChoice begin
    AC_NONE = 1
    AC_DT = 2
    AC_E0 = 3
    AC_DXT = 4
end

@enum JumpOffRate begin
    JR_FITTED = 1
    JR_ACTUAL = 2
end


@enum CalculationMode begin
    CM_DEMOGRAPHY = 1
    CM_JULIA = 2
end

@enum DataDimension begin
    DD_AGE = 1
    DD_YEAR = 2
end


@enum MortalityDataSource begin
    MDS_OBSERVED = 1
    MDS_FITTED = 2
    MDS_PREDICTED = 3
    MDS_CALCULATED = 4
end

@enum MortalityDataCategory begin
    MDC_EXPOSURES = 1
    MDC_DEATHS = 2
    MDC_APPROXIMATE_DEATHS = 3
    MDC_RATES = 4
    MDC_LOGRATES = 5
    MDC_LIFE_EXPECTANCIES = 6
    MDC_PARAMETERS = 7
    MDC_OTHER = 8
end

@enum UncertaintyMode begin
    UM_INNOVATION_ONLY = 1
    UM_INNOVDRIFT = 2
end





function sexmatch(sex::Sex, outcomes::Vararg{Union{T,Nothing},3}) where {T}
    outcome = @match sex begin
        $SEX_FEMALE => outcomes[2]
        $SEX_MALE => outcomes[3]
        _ => outcomes[1]
    end

    return outcome
end

function mdc_convert(s::Symbol)
    outcome = @match s begin
        :Exposures || :Ext => 1
        :Deaths || :Dxt => 2
        :ApproximateDeaths || :D̃xt => 3
        :Rates || :mxt => 4
        :LogRates || :logmxt => 5
        :LifeExpectancies || :ext => 6
        :Parameters => 7
        _ => 8
    end

    return MortalityDataCategory(outcome)

end

function mdc_label(input::Union{MortalityDataCategory,Symbol}, source::MortalityDataSource)
    datatype = input isa MortalityDataCategory ? input : mdc_convert(input)
    raw_label = @match Int(datatype) begin
        1 => "Central Exposed To Risk [Eᶜ(x,t)]"
        2 => "Observed Death Counts (calculated by HMD) [D(x,t)]"
        3 => "Approximated Death Counts (calculated by `demography` package as Eᶜ(x,t) × m(x,t)) [D̂(x,t)]"
        4 => "Central Rate Of Mortality [m(x,t)]"
        5 => "Natural Log Of Central Rate Of Mortality [ln(m(x,t))]"
        6 => "Mean Residual Lifetime (calculated by HMD) [e(x,t)]"
        7 => "Age×Year Parameter Values"
        8 => "Other"
    end


    if source == MDS_FITTED
        raw_label = "Fitted $raw_label"
    elseif source == MDS_PREDICTED
        raw_label = "Predicted $raw_label"
    elseif source == MDS_CALCULATED
        raw_label = "Caluclated $raw_label"
    end

    return raw_label

end


function mdc_shortlabel(input::Union{MortalityDataCategory,Symbol}, source::MortalityDataSource)
    datatype = input isa MortalityDataCategory ? input : md_convert(input, source)
    raw_label = @match Int(datatype) begin
        1 => "Eᶜ(x,t)"
        2 => "D(x,t)"
        3 => "D̃(x,t)"
        4 => "m(x,t)"
        5 => "ln(m(x,t))"
        6 => "e(x,t)"
        7 => "Age×Year Parameter Values"
        8 => "Other"
    end

    if source == MDS_FITTED
        raw_label = "Fitted $raw_label"
    elseif source == MDS_PREDICTED
        raw_label = "Predicted $raw_label"
    elseif source == MDS_CALCULATED
        raw_label = "Calculated $raw_label"
    end

    return raw_label
end


function a0_config(mode::CalculationMode, sex::Sex)
    if sex == SEX_BOTH
        throw(ArgumentError("Can't calculate a0 for both sexes"))
    end

    hmd_lb_m = [0.08307, 0.0230, 0]
    hmd_slopes_m = [0, 3.26021, -1.99545]
    hmd_intercepts_m = [0.29915, 0.02832, 0.14929]
    hmd_lb_f = [0.06891, 0.01724, 0]
    hmd_slopes_f = [0, 3.88089, -2.05527]
    hmd_intercepts_f = [0.01724, 0.04667, 0.14903]
    hmd = [(lb=hmd_lb_f, slope=hmd_slopes_f, intercept=hmd_intercepts_f), (lb=hmd_lb_m, slope=hmd_slopes_m, intercept=hmd_intercepts_m)]

    dem_lb_m = [0.107, 0]
    dem_lb_f = [0.107, 0]
    dem_slopes_m = [0, 2.684]
    dem_slopes_f = [0, 2.8]
    dem_intercepts_m = [0.33, 0.045]
    dem_intercepts_f = [0.35, 0.053]

    dem = [(lb=dem_lb_f, slope=dem_slopes_f, intercept=dem_intercepts_f), (lb=dem_lb_m, slope=dem_slopes_m, intercept=dem_intercepts_m)]

    results = [dem, hmd]

    return results[Int(mode)][Int(sex)-1]
end


# function Base.*(dc1::MortalityDataCategory, dc2::MortalityDataCategory)
#     v1 = Int(dc1)
#     v2 = Int(dc2)



# end