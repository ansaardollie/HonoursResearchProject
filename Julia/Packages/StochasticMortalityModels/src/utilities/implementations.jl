export ModelImplementation, bms_j, lc_j, lm_j, bms_r, lc_r, lm_r, bms, lc, lm
struct ModelImplementation
    adjustment::AdjustmentChoice
    jumpoff::JumpoffChoice
    calculation::CalculationChoice
end


const bms = ModelImplementation(AC_DXT, JR_FITTED, CC_JULIA)
const lc = ModelImplementation(AC_DT, JR_FITTED, CC_JULIA)
const lm = ModelImplementation(AC_E0, JR_ACTUAL, CC_JULIA)


const bms_j = ModelImplementation(AC_DXT, JR_FITTED, CC_JULIA)
const lc_j = ModelImplementation(AC_DT, JR_FITTED, CC_JULIA)
const lm_j = ModelImplementation(AC_E0, JR_ACTUAL, CC_JULIA)


const bms_r = ModelImplementation(AC_DXT, JR_FITTED, CC_DEMOGRAPHY)
const lc_r = ModelImplementation(AC_DT, JR_FITTED, CC_DEMOGRAPHY)
const lm_r = ModelImplementation(AC_E0, JR_ACTUAL, CC_DEMOGRAPHY)
