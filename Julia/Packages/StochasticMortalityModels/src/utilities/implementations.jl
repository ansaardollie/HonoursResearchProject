export ModelImplementation, bms, lc, lm
struct ModelImplementation
    adjustment::AdjustmentChoice
    jumpoff::JumpOffRate
end


const bms = ModelImplementation(AC_DXT, JR_FITTED)
const lc = ModelImplementation(AC_DT, JR_FITTED)
const lm = ModelImplementation(AC_E0, JR_ACTUAL)

ModelImplementation(AC_DXT, JR_ACTUAL)