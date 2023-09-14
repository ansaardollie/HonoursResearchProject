export AgePeriodData

struct AgePeriodData{T<:Real} <: AbstractMatrix{T}
    category::MortalityDataCategory
    source::MortalityDataSource
    label::String
    values::Matrix{T}
    ages::DataRange
    periods::DataRange
    rounding::Int
    summarised::Bool
    summaryfor::Optional{DataDimension}
end

@inline function AgePeriodData(
    category::Union{Symbol,MortalityDataCategory,Nothing},
    data::Matrix{T},
    ages::DataRange,
    periods::DataRange;
    rounding::Optional{Integer}=nothing,
    source::MortalityDataSource=MDS_OBSERVED,
    label::Optional{String}=nothing,
    summaryfor::Optional{DataDimension}=nothing
) where {T<:Real}

    if isnothing(summaryfor)
        length(ages) == size(data, 1) || throw(ArgumentError("Length of ages vector != # of rows in data"))
        length(periods) == size(data, 2) || throw(ArgumentError("Length of periods vector != # of columns in data"))
    elseif summaryfor == DD_AGE
        length(periods) == size(data, 2) || throw(ArgumentError("Length of periods vector != # of columns in data"))
    elseif summaryfor == DD_YEAR
        length(ages) == size(data, 1) || throw(ArgumentError("Length of ages vector != # of rows in data"))
    end

    label = isnothing(label) ?
            mdc_label(category, source) :
            label

    mdc_enum = category isa MortalityDataCategory ?
               category :
               category isa Symbol ? mdc_convert(category) :
               MDC_OTHER

    rounding = isnothing(rounding) ?
               [1, 3, 3, 6, 6, 2, 3, 6][Int(mdc_enum)] :
               rounding
    summarised = !isnothing(summaryfor)

    return AgePeriodData(
        mdc_enum,
        source,
        label,
        data,
        ages,
        periods,
        rounding,
        summarised,
        summaryfor
    )
end

@inline function AgePeriodData(
    category::Union{Symbol,MortalityDataCategory,Nothing},
    data::Matrix{T},
    ages::AbstractVector{Int},
    periods::AbstractVector{Int};
    rounding::Union{Int,Nothing}=nothing,
    source::MortalityDataSource=MDS_OBSERVED,
    label::Optional{String}=nothing,
    summaryfor::Optional{DataDimension}=nothing
) where {T<:Real}
    ar = agerange(ages)
    yr = yearrange(periods)

    return AgePeriodData(category, data, ar, yr; rounding=rounding, label=label, source=source, summaryfor=summaryfor)
end



function periods(atd::AgePeriodData{T}) where {T<:Real}
    return atd.periods
end

function ages(atd::AgePeriodData{T}) where {T<:Real}
    return atd.ages
end



function index_age(atd::AgePeriodData{T}, ages...) where {T<:Real}
    idx = map(a -> atd.ages.indexer[a], ages)
    return collect(idx)
end

function index_year(atd::AgePeriodData{T}, periods...) where {T<:Real}
    idx = map(y -> atd.periods.indexer[y], periods)
    return collect(idx)
end

function Base.axes(atd::AgePeriodData{T}) where {T<:Real}
    if isnothing(atd.summaryfor)
        return (atd.ages, atd.periods)
    else
        axis = atd.summaryfor == DD_AGE ? (atd.periods,) : (atd.ages,)
        return axis
    end
end


function Base.size(atd::AgePeriodData{T}) where {T<:Real}
    if !(atd.summarised)
        return (length(atd.ages), length(atd.periods))
    else
        # sizes = atd.summaryfor == DD_AGE ? (1, length(atd.periods)) : (length(atd.ages), 1)
        # return sizes
        return (length(atd.summaryfor == DD_AGE ? atd.periods : atd.ages),)
    end
end

function Base.IndexStyle(apd::AgePeriodData{T}) where {T<:Real}
    if apd.summarised
        return IndexLinear()
    else
        return IndexCartesian()
    end
end

function Base.LinearIndices(apd::AgePeriodData{T}) where {T<:Real}
    if apd.summarised
        return (apd.summaryfor == DD_AGE ? apd.periods.values : apd.ages.values)
    else
        return reshape(collect(Base.oneto(length(apd))), size(apd)...)
    end
end

@inline function Base.getindex(atd::AgePeriodData{T}, i::Int) where {T<:Real}
    @boundscheck begin
        if atd.summarised
            range = atd.summaryfor == DD_AGE ? atd.periods : atd.ages
            i in range || throw(BoundsError(range, i))
        else
            i <= length(atd) || throw(BoundsError(atd.values, i))
        end
    end

    if atd.summarised

        idx = atd.summaryfor == DD_AGE ? (1, atd.periods.indexer[i]) : (atd.ages.indexer[i], 1)

        return atd.values[idx...]

    else
        return atd.values[i]
    end
end

function Base.length(atd::AgePeriodData{T}) where {T}
    return length(atd.values)
end
function Base.IteratorSize(atd::AgePeriodData{T}) where {T}
    if atd.summarised
        return Base.HasLength()
    else
        return Base.HasShape{2}()
    end
end

@inline function Base.getindex(atd::AgePeriodData{T}, age::Int, year::Int) where {T<:Real}
    @boundscheck begin
        (haskey(atd.ages.indexer, age) && haskey(atd.periods.indexer, year)) || throw(BoundsError(atd, (age=age, year=year)))
    end

    @inbounds row = atd.ages.indexer[age]
    @inbounds col = atd.periods.indexer[year]
    return @inbounds atd.values[row, col]
end



@inline function Base.getindex(atd::AgePeriodData{T}, I...) where {T<:Real}
    if I[1] isa CartesianIndex
        idx_a = I[1].I[1]
        idx_y = I[1].I[2]
        return atd[idx_a, idx_y]
    end



    @boundscheck begin
        if !atd.summarised
            length(I) <= 2 || throw(BoundsError(atd, I))
            I[1] isa Colon || issubset(I[1], atd.ages.values) || throw(BoundsError(atd, I[1]))
            I[2] isa Colon || issubset(I[2], atd.periods.values) || throw(BoundsError(atd, I[2]))
        else
            length(I) == 1 || throw(BoundsError("Can't index summarised AgePeriodData with more than one index."))
            range = atd.summaryfor == DD_AGE ? atd.periods : atd.ages
            I[1] isa Colon || issubset(I[1], range) || throw(BoundsError(atd, I[1]))
        end
    end

    if !atd.summarised

        idx_a = I[1]
        idx_y = I[2]

        row_idx = idx_a isa Colon ? idx_a : index_age(atd, idx_a...)
        col_idx = idx_y isa Colon ? idx_y : index_year(atd, idx_y...)


        # (idx_a isa Colon) || (size(idx_a) == () && (row_idx = [row_idx]; true))
        # (idx_y isa Colon) || (size(idx_y) == () && (col_idx = [col_idx]; true))

        return @inbounds AgePeriodData(
            atd.category,
            atd.values[row_idx, col_idx],
            atd.ages.values[row_idx],
            atd.periods.values[col_idx];
            rounding=atd.rounding,
            source=atd.source,
            summaryfor=atd.summaryfor
        )
    elseif atd.summaryfor == DD_AGE

        idx_y = I[1]

        row_idx = [1]
        col_idx = idx_y isa Colon ? idx_y : index_year(atd, idx_y...)

        return @inbounds AgePeriodData(
            atd.category,
            atd.values[row_idx, col_idx],
            atd.ages,
            atd.periods.values[col_idx];
            rounding=atd.rounding,
            source=atd.source,
            summaryfor=atd.summaryfor
        )

    elseif atd.summaryfor == DD_YEAR
        idx_a = I[1]

        row_idx = idx_a isa Colon ? idx_a : index_age(atd, idx_a...)
        col_idx = [1]

        return @inbounds AgePeriodData(
            atd.category,
            atd.values[row_idx, col_idx],
            atd.ages.values[row_idx],
            atd.periods;
            rounding=atd.rounding,
            source=atd.source,
            summaryfor=atd.summaryfor
        )

    end
end

@inline function Base.setindex!(atd::AgePeriodData{T}, V::T, i::Int) where {T<:Real}
    @boundscheck begin
        if atd.summarised
            range = atd.summaryfor == DD_AGE ? atd.periods : atd.ages
            i in range || throw(BoundsError(range, i))
        else
            i <= length(atd) || throw(BoundsError(atd.values, i))
        end
    end

    if atd.summarised

        @inbounds idx = atd.summaryfor == DD_AGE ? (1, atd.periods.indexer[i]) : (atd.ages.indexer[i], 1)

        @inbounds atd.values[idx...] = V
    else
        @inbounds atd.values[i] = V
    end
    return atd
end

@inline function Base.setindex!(atd::AgePeriodData{T}, V::T, age::Int, year::Int) where {T<:Real}
    @boundscheck begin
        (haskey(atd.ages.indexer, age) && haskey(atd.periods.indexer, year)) || throw(BoundsError(atd, (age=age, year=year)))
    end

    @inbounds row = atd.ages.indexer[age]
    @inbounds col = atd.periods.indexer[year]

    @inbounds atd.values[row, col] = V

    return atd

end

@inline function Base.setindex!(atd::AgePeriodData{T}, V::AbstractArray{T}, I...) where {T<:Real}

    if I[1] isa CartesianIndex
        idx_a = I[1].I[1]
        idx_y = I[1].I[2]
        return atd[idx_a, idx_y] = v
    end

    @boundscheck begin
        if !atd.summarised
            length(I) <= 2 || throw(BoundsError(atd, I))
            I[1] isa Colon || issubset(I[1], atd.ages.values) || throw(BoundsError(atd, I[1]))
            I[2] isa Colon || issubset(I[2], atd.periods.values) || throw(BoundsError(atd, I[2]))
        else
            length(I) == 1 || throw(BoundsError("Can't index summarised AgePeriodData with more than one index."))
            range = atd.summaryfor == DD_AGE ? atd.periods : atd.ages
            I[1] isa Colon || issubset(I[1], range) || throw(BoundsError(atd, I[1]))
        end
    end

    if !atd.summarised
        idx_a = I[1]
        idx_y = I[2]

        row_idx = idx_a isa Colon ? idx_a : index_age(atd, idx_a...)
        col_idx = idx_y isa Colon ? idx_y : index_year(atd, idx_y...)


        atd.values[row_idx, col_idx] = V
    elseif atd.summaryfor == DD_AGE
        idx_a = 1
        idx_y = I[1]
        row_idx = 1
        col_idx = idx_y isa Colon ? idx_y : index_year(atd, idx_y...)

        atd.values[row_idx, col_idx] = V

    elseif atd.summaryfor == DD_YEAR
        idx_a = I[1]
        idx_y = 1
        row_idx = idx_a isa Colon ? idx_a : index_age(atd, idx_a...)
        col_idx = 1

        atd.values[row_idx, col_idx] = V
    end



    return atd

end

function Base.iterate(atd::AgePeriodData{T}) where {T<:Real}
    if atd.summarised
        range = atd.summaryfor == DD_AGE ? atd.periods : atd.ages

        out = iterate(range)
        if isnothing(out)
            return nothing
        else
            return (atd[out[1]], out[2])
        end
    else
        years = atd.periods
        ages = atd.ages
        ia = iterate(ages)
        iy = iterate(years)
        if length(years) == 0 && length(ages) == 0
            return nothing
        else
            age_idx = ia[1]
            year_idx = iy[1]
            age_state = ia[2]
            year_state = iy[2]
            return (atd[age_idx, year_idx], ((age_state, year_state), (age_idx, year_idx)))
        end
    end

end


function Base.iterate(atd::AgePeriodData{T}, state::Int) where {T<:Real}
    if !atd.summarised
        return nothing
    end

    range = atd.summaryfor == DD_AGE ? atd.periods : atd.ages

    out = iterate(range, state)
    if isnothing(out)
        return nothing
    else
        return (atd[out[1]], out[2])
    end
end


function Base.iterate(atd::AgePeriodData{T}, state::Tuple{Tuple{Int,Int},Tuple{Int,Int}}) where {T<:Real}
    if atd.summarised
        return nothing
    end

    years = atd.periods
    ages = atd.ages
    ia = iterate(ages, state[1][1])
    iy = iterate(years, state[1][2])
    if isnothing(ia) && isnothing(iy)
        return nothing
    elseif isnothing(ia)
        ia = iterate(ages)
        age_idx = ia[1]
        year_idx = iy[1]
        age_state = ia[2]
        year_state = iy[2]
        return (atd[age_idx, year_idx], ((age_state, year_state), (age_idx, year_idx)))
    else
        age_idx = ia[1]
        year_idx = state[2][2]
        age_state = ia[2]
        year_state = state[1][2]
        return (atd[age_idx, year_idx], ((age_state, year_state), (age_idx, year_idx)))
    end
end


function Base.firstindex(atd::AgePeriodData{T}) where {T}
    if atd.summarised
        range = atd.summaryfor == DD_AGE ? atd.periods : atd.ages
        return Base.firstindex(range)
    else
        age_idx = Base.firstindex(atd.ages)
        year_idx = Base.firstindex(atd.periods)

        return CartesianIndex(age_idx, year_idx)
    end
end

function Base.lastindex(atd::AgePeriodData{T}) where {T}
    if atd.summarised
        range = atd.summaryfor == DD_AGE ? atd.periods : atd.ages
        return Base.lastindex(range)
    else
        age_idx = Base.lastindex(atd.ages)
        year_idx = Base.lastindex(atd.periods)

        return CartesianIndex(age_idx, year_idx)
    end
end

function generate_formatter(data::Matrix{T}, defaultrounding::Union{Int,Nothing}=nothing) where {T<:Real}



    non_zeros = filter(x -> x != zero(T), data)
    if isempty(non_zeros)
        return (formatter=ft_round(0), arm=1, max=0, min=0, magrange=0)
    end

    max_exponent = Int(ceil(log10(maximum(abs.(non_zeros)))))
    min_exponent = Int(floor(log10(minimum(abs.(non_zeros)))))

    max_value = maximum(abs.(data))
    mv_exponent = Int(ceil(log10(max_value)))


    range_magnitude = max_exponent - min_exponent

    arm = 0
    if range_magnitude <= 5 || !(isnothing(defaultrounding))
        arm = 1
    elseif -10 < min_exponent <= max_exponent < 10
        arm = 2
    elseif range_magnitude >= 20
        arm = 3
    elseif abs(min_exponent) < max_exponent
        arm = 4
    else
        arm = 5
    end

    function format_func(value::T, row::Number, col::Number)
        if value == zero(T)
            return "0"
        end
        left_digits = abs(value) >= 1 ? Int(ceil(log10(abs(value)))) : 1
        exponent = Int(floor(log10(abs(value))))
        sexp = exponent < 0 ? "$exponent" : "+$exponent"
        mantissa = value / 10.0^(exponent)
        small_rep = value / 10.0^(max_exponent)
        small_rep_ld = abs(small_rep) >= 1 ? Int(ceil(log10(abs(small_rep)))) : 1
        large_rep = value / 10.0^(min_exponent)
        large_rep_ld = abs(large_rep) >= 1 ? Int(ceil(log10(abs(large_rep)))) : 1

        outcome = @match arm begin
            1 => round(value; digits=defaultrounding)
            2 => round(value; digits=9 - left_digits)
            3 => "$(round(mantissa; digits=6))e$sexp"
            4 => round(small_rep; digits=9 - small_rep_ld)
            5 => round(large_rep; digits=9 - large_rep_ld)
            _ => round(value; digit=9 - left_digits)
        end

        if arm != 3
            decimals = @match arm begin
                1 => defaultrounding
                2 => 9 - left_digits
                4 => 9 - small_rep_ld
                5 => 9 - large_rep_ld
                _ => 9 - left_digits
            end
            fmtString = "%.$(decimals)f"
            outcome = Printf.format(Printf.Format(fmtString), outcome)
        end

        parts = split(outcome, letter -> letter == '.' || letter == 'e')

        lhs = reverse(parts[1])
        rhs = parts[2]

        lseperated = reverse(join(
            [(i % 3 == 0) ? "$(lhs[i]) " : "$(lhs[i])" for i in eachindex(lhs)]
        ))
        rseperated = join([(i % 3 == 0) ? "$(rhs[i]) " : "$(rhs[i])" for i in eachindex(rhs)])

        outcome = strip("$(lseperated).$(rseperated)")
        if length(parts) >= 3
            outcome = "$(outcome)e$(parts[3])"
        end

        return outcome
    end

    return (formatter=format_func, arm=arm, max=max_exponent, min=min_exponent, magrange=range_magnitude)
end


function generate_summary(label::Union{String,Nothing}, ages::DataRange, periods::DataRange, fmt_results::Any)
    titlebuf = IOBuffer()
    println(titlebuf)
    println(titlebuf)
    if !isnothing(label)
        println(titlebuf, "Dataset: $(isnothing(label) ? "Unknown" : label)")
    end

    println(titlebuf, "Ages: $(ages)")
    println(titlebuf, "Years: $(periods)")

    if fmt_results.arm >= 4
        exp = fmt_results.arm == 4 ? fmt_results.max : fmt_results.min
        exp = (exp > 0) ? "+$exp" : "$exp"
        println(titlebuf)
        println(titlebuf, "Values in units of 10^($(exp))")
    end

    title_out = String(take!(titlebuf))

    return title_out
end


function Base.show(io::IO, ::MIME"text/plain", atd::AgePeriodData{T}) where {T<:Real}
    periods = atd.periods.values
    ages = atd.ages.values


    (ff, format_results...) = generate_formatter(atd.values, atd.rounding)


    summary = generate_summary(atd.label, atd.ages, atd.periods, format_results)

    headers = atd.summaryfor == DD_YEAR ? [Symbol("$(atd.periods)")] : map(Symbol, periods)
    rows = atd.summaryfor == DD_AGE ? [Symbol("$(atd.ages)")] : ages

    apd_conf = table_config(io;
        headers=headers,
        rows=rows,
        row_label_title="""Age \\ Year""",
        formatters=ff,
        title=summary
    )


    pretty_table(io, atd.values; apd_conf...)
end



function Base.:+(atd1::AgePeriodData, atd2::AgePeriodData)
    size(atd1) == size(atd2) || throw(DimensionMismatch("Sizes not compatible"))

    sum(atd1.ages .== atd2.ages) == length(atd1.ages) || throw(DimensionMismatch("Ages $(atd1.ages) and $(atd2.ages) not compatible"))


    sum(atd1.periods .== atd2.periods) == length(atd1.periods) || throw(DimensionMismatch("Periods $(atd1.periods) and $(atd2.periods) not compatible"))

    result = atd1.values .+ atd2.values
    res_source = atd1.source == atd2.source ? atd1.source : MDS_CALCULATED
    res_round = atd1.category == atd2.category ? atd1.rounding : 6

    if (atd1.category == atd2.category) & (atd1.category != MDC_PARAMETERS)
        return AgePeriodData(atd1.category, result, atd1.ages, atd1.periods; source=res_source)
    elseif atd1.category == atd2.category
        n = "$(atd1.label) + $(atd2.label)"
        return AgePeriodData(MDC_PARAMETERS, result, atd1.ages, atd2.periods; label=n, source=res_source)
    else
        n = "$(mds_shortlabel(atd1.category, atd1.source)) + $(mds_shortlabel(atd2.category,atd2.source))"
        return AgePeriodData(MDC_OTHER, result, atd1.ages, atd2.periods; label=n, source=res_source, rounding=res_round)
    end
end

function Base.:*(atd1::AgePeriodData, atd2::AgePeriodData)
    size(atd1) == size(atd2) || throw(DimensionMismatch("Sizes not compatible"))

    sum(atd1.ages .== atd2.ages) == length(atd1.ages) || throw(DimensionMismatch("Ages $(atd1.ages) and $(atd2.ages) not compatible"))


    sum(atd1.periods .== atd2.periods) == length(atd1.periods) || throw(DimensionMismatch("Periods $(atd1.periods) and $(atd2.periods) not compatible"))

    result = atd1.values .* atd2.values
    res_source = atd1.source == atd2.source ? atd1.source : MDS_CALCULATED
    res_round = atd1.category == atd2.category ? atd1.rounding : 6

    if (atd1.category == atd2.category) & (atd1.category != MDC_PARAMETERS)
        return AgePeriodData(atd1.category, result, atd1.ages, atd1.periods; source=res_source)
    elseif atd1.category == atd2.category
        n = "$(atd1.label) × $(atd2.label)"
        return AgePeriodData(MDC_PARAMETERS, result, atd1.ages, atd2.periods; label=n, source=res_source)
    else
        n = "$(mdc_shortlabel(atd1.category, atd1.source)) × $(mdc_shortlabel(atd2.category,atd2.source))"
        return AgePeriodData(MDC_OTHER, result, atd1.ages, atd2.periods; label=n, source=res_source, rounding=res_round)
    end
end

# MyDims = Union{Dims,Tuple{Vararg{DataRange,N}}} where {N}
MyDims = Tuple{Vararg{Union{DataRange,Int},N}} where {N}

function Base.similar(atd::AgePeriodData{T}, ::Type{S}, dim::MyDims) where {T<:Real,S<:Real}

    if atd.summarised

        # ages = atd.summaryfor == DD_AGE ? atd.ages : dim[1]
        if atd.summaryfor == DD_AGE
            ages = atd.ages

            year_ranges = filter(d -> d isa DataRange && d.type == DD_YEAR, dim)
            potential_ranges = [atd.periods, year_ranges...]
            range_lengths = map(pr -> length(pr), potential_ranges)
            (max_length, idx) = findmax(range_lengths)
            range = potential_ranges[idx]
            periods = range
        else
            years = atd.periods

            age_ranges = filter(d -> d isa DataRange && d.type == DD_AGE, dim)
            potential_ranges = [atd.ages, age_ranges...]
            range_lengths = map(pr -> length(pr), potential_ranges)
            (max_length, idx) = findmax(range_lengths)
            range = potential_ranges[idx]
            ages = range
        end
        n = length(range)
        size = atd.summaryfor == DD_AGE ? (1, n) : (n, 1)

        res = Matrix{S}(undef, size...)
        mdc_enum = MDC_OTHER
        src = atd.source
        rounding = atd.rounding
        n = "f($(atd.category == MDC_PARAMETERS ? atd.label : mdc_shortlabel(atd.category,atd.source)))"
        return AgePeriodData(mdc_enum, res, ages, periods; source=src, rounding=rounding, label=n, summaryfor=atd.summaryfor)

    elseif dim[1] isa Int
        ages = atd.ages
        periods = dim[2]
        n = length(periods)
        size = (1, n)
        res = Matrix{S}(undef, size...)
        mdc_enum = MDC_OTHER
        src = atd.source
        rounding = atd.rounding
        n = "f($(atd.category == MDC_PARAMETERS ? atd.label : mdc_shortlabel(atd.category,atd.source)))"
        return AgePeriodData(mdc_enum, res, ages, periods; source=src, rounding=rounding, label=n, summaryfor=DD_AGE)
    elseif dim[2] isa Int
        ages = dim[1]
        periods = atd.periods
        n = length(ages)
        size = (n, 1)
        res = Matrix{S}(undef, size...)
        mdc_enum = MDC_OTHER
        src = atd.source
        rounding = atd.rounding
        n = "f($(atd.category == MDC_PARAMETERS ? atd.label : mdc_shortlabel(atd.category,atd.source)))"
        return AgePeriodData(mdc_enum, res, ages, periods; source=src, rounding=rounding, label=n, summaryfor=DD_YEAR)
    else
        ages = dim[1]
        n = length(ages)
        periods = dim[2]
        m = length(periods)

        res = Matrix{S}(undef, n, m)
        mdc_enum = MDC_OTHER
        src = atd.source
        rounding = atd.rounding
        n = "f($(atd.category == MDC_PARAMETERS ? atd.label : mdc_shortlabel(atd.category,atd.source)))"
        return AgePeriodData(mdc_enum, res, ages, periods; source=src, rounding=rounding, label=n)
    end

end


function Base.reduced_indices(apd::AgePeriodData{T}, region::UnitRange{Int64}) where {T}
    if length(apd.periods) == 1
        return (1, apd.periods)
    elseif length(apd.ages) == 1
        return (apd.ages, 1)
    end
end