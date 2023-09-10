# include("./enums.jl")
# include("./datarange.jl")

# using PrettyTables
# using Printf
# struct AgeYearData{T<:Real} <: AbstractMatrix{T}
#     rawtype::Union{ObservedMortalityDataType,PredictedMortalityDataType,FittedMortalityDataType,Nothing}
#     label::Union{String,Nothing}
#     data::Matrix{T}
#     years::DataRange
#     ages::DataRange
#     year_indexer::Dict{Int,Int}
#     age_indexer::Dict{Int,Int}
#     default_rounding::Int
#     datakind::MortalityDataType
#     @inline function AgeYearData(rt::Any, data::Matrix{T}, ages::DataRange, years::DataRange; rounding::Union{Integer,Nothing}=nothing, datakind::MortalityData=MDT_OBSERVED) where {T<:Real}

#         length(ages) == size(data, 1) || throw(ArgumentError("Length of ages vector != # of rows in data"))
#         length(years) == size(data, 2) || throw(ArgumentError("Length of years vector != # of columns in data"))
#         yid = Dict(zip(years, eachindex(years)))
#         aid = Dict(zip(ages, eachindex(ages)))

#         label = rt isa String ? rt : md_label(rt, datakind)
#         rawtype_enum = rt isa String ?
#                        [OMD_PARAMETERS, FMD_PARAMETERS, PMD_PARAMETERS][Int(datakind)] :
#                        rt isa MortalityData ?
#                        rt :
#                        md_convert(rt, datakind)

#         df_rounding = isnothing(rounding) ? [1, 3, 3, 6, 6, 2, 6][Int(rt)] : rounding

#         return new{T}(rawtype_enum, label, data, years, ages, yid, aid, df_rounding, datakind)
#     end
# end

# function AgeYearData(
#     rawtype::Union{Nothing,String,Symbol,ObservedMortalityDataType,FittedMortalityDataType,PredictedMortalityDataType},
#     data::Matrix{T},
#     ages::DataRange,
#     years::DataRange;
#     rounding::Union{Int,Nothing}=nothing,
#     datakind::MortalityData=MDT_OBSERVED
# ) where {T<:Real}

#     println("Inside CTOR1")
#     # @boundscheck begin
#     length(ages) == size(data, 1) || throw(ArgumentError("Length of ages vector != # of rows in data"))
#     length(years) == size(data, 2) || throw(ArgumentError("Length of years vector != # of columns in data"))
#     # end
#     yid = Dict(zip(years, eachindex(years)))
#     aid = Dict(zip(ages, eachindex(ages)))


#     label = rawtype isa String ? rawtype : md_label(rawtype, datakind)
#     rawtype_enum = rawtype isa String ?
#                    [OMD_PARAMETERS, FMD_PARAMETERS, PMD_PARAMETERS][Int(datakind)] :
#                    rawtype isa MortalityData ?
#                    rawtype :
#                    md_convert(rawtype, datakind)

#     df_rounding = isnothing(rounding) ? [1, 3, 3, 6, 6, 2, 6][Int(rawtype)] : rounding

#     return AgeYearData{T}(rawtype_enum, label, data, years, ages, yid, aid, df_rounding, datakind)
# end
# @inline function AgeYearData(
#     rawtype::Union{Nothing,String,Symbol,ObservedMortalityDataType,FittedMortalityDataType,PredictedMortalityDataType},
#     data::Matrix{T},
#     ages::Vector{Int},
#     years::Vector{Int};
#     rounding::Union{Int,Nothing}=nothing,
#     datakind::MortalityData=MDT_OBSERVED
# ) where {T<:Real}

#     println("Inside CTOR2")
#     ydr = DataRange(years, DD_YEAR)
#     adr = DataRange(ages, DD_AGE)

#     return AgeYearData(rawtype, data, adr, ydr; rounding=rounding, datakind=datakind)
# end
# @inline function AgeYearData(
#     data::Matrix{T},
#     ages::AbstractVector{Int},
#     years::AbstractVector{Int};
#     rounding::Union{Int,Nothing}=nothing,
#     datakind::MortalityData=MDT_OBSERVED
# ) where {T<:Real}
#     println("Inside CTOR3")
#     return AgeYearData(nothing, data, ages, years; rounding=rounding, datakind=datakind)
# end


# function generate_formatter(data::Matrix{T}, defaultrounding::Union{Int,Nothing}=nothing) where {T<:Real}

#     non_zeros = filter(x -> x != zero(T), data)
#     max_exponent = Int(ceil(log10(maximum(abs.(non_zeros)))))
#     min_exponent = Int(floor(log10(minimum(abs.(non_zeros)))))

#     max_value = maximum(abs.(data))
#     mv_exponent = Int(ceil(log10(max_value)))
#     # mv_round = isnothing(defaultrounding) ? 9 - mv_exponent : defaultrounding

#     range_magnitude = max_exponent - min_exponent

#     arm = 0
#     if range_magnitude <= 5 || !(isnothing(defaultrounding))
#         arm = 1
#     elseif -10 < min_exponent <= max_exponent < 10
#         arm = 2
#     elseif range_magnitude >= 20
#         arm = 3
#     elseif abs(min_exponent) < max_exponent
#         arm = 4
#     else
#         arm = 5
#     end

#     function format_func(value::T, row::Number, col::Number)
#         if value == zero(T)
#             return "0"
#         end
#         left_digits = abs(value) >= 1 ? Int(ceil(log10(abs(value)))) : 1
#         exponent = Int(floor(log10(abs(value))))
#         sexp = exponent < 0 ? "$exponent" : "+$exponent"
#         mantissa = value / 10.0^(exponent)
#         small_rep = value / 10.0^(max_exponent)
#         small_rep_ld = abs(small_rep) >= 1 ? Int(ceil(log10(abs(small_rep)))) : 1
#         large_rep = value / 10.0^(min_exponent)
#         large_rep_ld = abs(large_rep) >= 1 ? Int(ceil(log10(abs(large_rep)))) : 1

#         outcome = @match arm begin
#             1 => round(value; digits=defaultrounding)
#             2 => round(value; digits=9 - left_digits)
#             3 => "$(round(mantissa; digits=6))e$sexp"
#             4 => round(small_rep; digits=9 - small_rep_ld)
#             5 => round(large_rep; digits=9 - large_rep_ld)
#             _ => round(value; digit=9 - left_digits)
#         end

#         if arm != 3
#             decimals = @match arm begin
#                 1 => defaultrounding
#                 2 => 9 - left_digits
#                 4 => 9 - small_rep_ld
#                 5 => 9 - large_rep_ld
#                 _ => 9 - left_digits
#             end
#             fmtString = "%.$(decimals)f"
#             outcome = Printf.format(Printf.Format(fmtString), outcome)
#         end

#         parts = split(outcome, letter -> letter == '.' || letter == 'e')

#         lhs = reverse(parts[1])
#         rhs = parts[2]

#         lseperated = reverse(join(
#             [(i % 3 == 0) ? "$(lhs[i]) " : "$(lhs[i])" for i in eachindex(lhs)]
#         ))
#         rseperated = join([(i % 3 == 0) ? "$(rhs[i]) " : "$(rhs[i])" for i in eachindex(rhs)])

#         outcome = strip("$(lseperated).$(rseperated)")
#         if length(parts) >= 3
#             outcome = "$(outcome)e$(parts[3])"
#         end

#         return outcome
#     end

#     return (formatter=format_func, arm=arm, max=max_exponent, min=min_exponent, magrange=range_magnitude)
# end

# function generate_summary(label::Union{String,Nothing}, ages::DataRange, years::DataRange, fmt_results::Any)
#     titlebuf = IOBuffer()
#     println(titlebuf)
#     println(titlebuf)
#     if !isnothing(label)
#         println(titlebuf, "Dataset: $(isnothing(label) ? "Unknown" : label)")
#     end

#     println(titlebuf, "$(ages)")
#     println(titlebuf, "$(years)")

#     if fmt_results.arm >= 4
#         exp = fmt_results.arm == 4 ? fmt_results.max : fmt_results.min
#         exp = (exp > 0) ? "+$exp" : "$exp"
#         println(titlebuf)
#         println(titlebuf, "Values in units of 10^($(exp))")
#     end

#     title_out = String(take!(titlebuf))

#     return title_out
# end

# function years(ayd::AgeYearData{T}) where {T<:Real}
#     return ayd.years
# end

# function ages(ayd::AgeYearData{T}) where {T<:Real}
#     return ayd.ages
# end

# function index_age(ayd::AgeYearData{T}, ages...) where {T<:Real}
#     idx = map(a -> ayd.ages[a], ages)
#     return idx
# end

# function index_year(ayd::AgeYearData{T}, years...) where {T<:Real}
#     idx = map(y -> ayd.year_indexer[y], years)
#     return idx
# end

# function Base.size(ayd::AgeYearData{T}) where {T<:Real}
#     return (length(ayd.ages), length(ayd.years))
# end

# function Base.IndexStyle(ayd::AgeYearData{T}) where {T<:Real}
#     return IndexCartesian()
# end

# function Base.axes(ayd::AgeYearData{T}) where {T<:Real}
#     return (ayd.ages, ayd.years)
# end


# @inline function Base.getindex(ayd::AgeYearData{T}, age::Int, year::Int) where {T<:Real}
#     @boundscheck begin
#         (haskey(ayd.age_indexer, age) && haskey(ayd.year_indexer, year)) || throw(BoundsError(ayd, (age=age, year=year)))
#     end

#     @inbounds row = ayd.age_indexer[age]
#     @inbounds col = ayd.year_indexer[year]
#     return @inbounds ayd.data[row, col]
# end


# @inline function Base.getindex(ayd::AgeYearData{T}, I...) where {T<:Real}
#     if I[1] isa CartesianIndex
#         idx_a = I[1].I[1]
#         idx_y = I[1].I[2]
#         return Base.getindex(ayd, idx_a, idx_y)
#     end

#     @boundscheck begin
#         length(I) <= 2 || throw(BoundsError(ayd, I))
#         idx_a = I[1]
#         idx_y = I[2]
#         idx_a isa Colon || issubset(idx_a, ayd.ages.values) || throw(BoundsError(ayd, I[1]))
#         idx_y isa Colon || issubset(idx_y, ayd.years.values) || throw(BoundsError(ayd, I[2]))
#     end

#     row_idx = idx_a isa Colon ? idx_a : map(age -> ayd.age_indexer[age], idx_a)
#     col_idx = idx_y isa Colon ? idx_y : map(year -> ayd.year_indexer[year], idx_y)


#     (idx_a isa Colon) || (size(idx_a) == () && (row_idx = [row_idx]; true))
#     (idx_y isa Colon) || (size(idx_y) == () && (col_idx = [col_idx]; true))

#     return @inbounds AgeYearData(
#         ayd.rawtype,
#         ayd.data[row_idx, col_idx],
#         ayd.ages.values[row_idx],
#         ayd.years.values[col_idx];
#         rounding=ayd.default_rounding
#     )
# end

# @inline function Base.setindex!(ayd::AgeYearData{T}, V::T, age::Int, year::Int) where {T<:Real}
#     @boundscheck begin
#         (haskey(ayd.age_indexer, age) && haskey(ayd.year_indexer, year)) || throw(BoundsError(ayd, (age=age, year=year)))
#     end

#     @inbounds row = ayd.age_indexer[age]
#     @inbounds col = ayd.year_indexer[year]

#     @inbounds ayd.data[row, col] = V

#     return @inbounds ayd.data[row, col]

# end

# @inline function Base.setindex!(ayd::AgeYearData{T}, V::AbstractArray{T}, I) where {T<:Real}

#     if I[1] isa CartesianIndex
#         idx_a = I[1].I[1]
#         idx_y = I[1].I[2]
#         return Base.setindex!(ayd, V[1], idx_a, idx_y)
#     end

#     @boundscheck begin
#         length(I) <= 2 || throw(BoundsError(ayd, I))
#         idx_a = I[1]
#         idx_y = I[2]
#         idx_a isa Colon || issubset(idx_a, ayd.ages.values) || throw(BoundsError(ayd, I[1]))
#         idx_y isa Colon || issubset(idx_y, ayd.years.values) || throw(BoundsError(ayd, I[2]))
#     end

#     row_idx = idx_a isa Colon ? idx_a : map(age -> ayd.age_indexer[age], idx_a)
#     col_idx = idx_y isa Colon ? idx_y : map(year -> ayd.year_indexer[year], idx_y)


#     (idx_a isa Colon) || (size(idx_a) == () && (row_idx = [row_idx]; true))
#     (idx_y isa Colon) || (size(idx_y) == () && (col_idx = [col_idx]; true))

#     return Base.setindex!(ayd.data, V, row_idx, col_idx)

# end



# function Base.show(io::IO, ::MIME"text/plain", ayd::AgeYearData{T}) where {T<:Real}
#     years = ayd.years.values
#     ages = ayd.ages.values


#     (ff, format_results...) = generate_formatter(ayd.data, ayd.default_rounding)


#     summary = generate_summary(ayd.label, ayd.ages, ayd.years, format_results)


#     rlh_crayon = crayon"fg:white bold italics"
#     title_crayon = crayon"bold"
#     rl_crayon = crayon"fg:white italics"
#     ch_crayon = crayon"fg:white italics"

#     conf = (
#         alignment=T <: Integer ? :c : :r,
#         backend=Val(:text),
#         formatters=ff,
#         header=map(Symbol, years),
#         header_alignment=:c,
#         limit_printing=false,
#         row_labels=ages,
#         row_label_alignment=:c,
#         row_label_column_title="""Age \\ Year""",
#         show_header=true,
#         show_row_number=false,
#         show_subheader=false,
#         title=summary,
#         vcrop_mode=:middle,
#         title_alignment=:l,
#         ellipsis_line_skip=0,
#         crop=:both,
#         linebreaks=true,
#         equal_columns_width=false,
#         title_same_width_as_table=true,
#         row_label_header_crayon=rlh_crayon,
#         title_crayon=title_crayon,
#         row_label_crayon=rl_crayon,
#         header_crayon=ch_crayon)

#     pretty_table(io, ayd.data; conf...)
# end


# function Base.:+(ayd1::AgeYearData, ayd2::AgeYearData)
#     size(ayd1) == size(ayd2) || throw(DimensionMismatch("Sizes not compatible"))

#     sum(ayd1.ages .== ayd2.ages) == length(ayd1.ages) || throw(DimensionMismatch("Ages $(ayd1.ages) and $(ayd2.ages) not compatible"))


#     sum(ayd1.years .== ayd2.years) == length(ayd1.years) || throw(DimensionMismatch("Years $(ayd1.years) and $(ayd2.years) not compatible"))

#     result = ayd1.data .+ ayd2.data

#     if (ayd1.rawtype == ayd2.rawtype) & (ayd1.rawtype != RD_PARAMETERS)
#         return AgeYearData(ayd1.rawtype, result, ayd1.ages, ayd1.years)
#     elseif ayd1.rawtype == ayd2.rawtype
#         n = "$(ayd1.label) + $(ayd2.label)"
#         return AgeYearData(n, result, ayd1.ages, ayd2.years)
#     else
#         n = "$(md_shortlabel(ayd1.rawtype, ayd1.datakind)) + $(md_shortlabel(ayd2.rawtype,ayd2.datakind))"
#         return AgeYearData(n, result, ayd1.ages, ayd2.years)
#     end
# end



# MyDims = Union{Dims,Tuple{Vararg{DataRange,N}}} where {N}

# function Base.similar(ps::AgeYearData{T}, el::Type{S}, dim::MyDims) where {T<:Real,S<:Real}

#     ages = dim[1]
#     n = length(ages)
#     years = dim[2]
#     m = length(years)

#     res = Matrix{S}(undef, n, m)
#     return AgeYearData("f($(ps.rawtype == RD_PARAMETERS ? ps.label : md_shortlabel(ps.rawtype,ps.datakind)))", res, ages, years; datakind=ps.datakind)

# end

# ay1 = AgeYearData(OMD_DEATHS, randn(4, 4), agerange(0:3), yearrange(0:3))