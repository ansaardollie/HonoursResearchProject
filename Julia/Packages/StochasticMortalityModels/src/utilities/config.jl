export table_config

ft_identity = (v, i, j) -> v


function table_config(io::IO;
    title=nothing, headers, rows=nothing, row_label_title=nothing,
    alignment=:r,
    backend=Val(:text),
    formatters=ft_identity,
    header_alignment=:c,
    limit_printing=false,
    row_label_alignment=:c,
    show_header=true,
    show_row_number=false,
    show_subheader=true,
    ellipsis_line_skip=0,
    linebreaks=true,
    crop=nothing,
    equal_columns_width=false,
    title_same_width_as_table=false,
    title_alignment=:l,
    hlines=nothing
)

    rlh_crayon = crayon"fg:white bold italics"
    title_crayon = crayon"bold"
    rl_crayon = crayon"fg:white italics"
    ch_crayon = crayon"fg:white italics"

    isfile = get(io, :file, false)::Bool


    conf = (
        title_alignment=title_alignment,
        backend=backend,
        formatters=formatters,
        limit_printing=limit_printing,
        alignment=alignment,
        linebreaks=linebreaks,
        header=headers,
        tf=isfile ? tf_compact : tf_unicode,
        show_subheader=show_subheader,
        header_alignment=header_alignment,
        row_label_alignment=row_label_alignment,
        show_header=show_header,
        show_row_number=show_row_number,
        crop=!isnothing(crop) ? crop : isfile ? :none : :both,
        vcrop_mode=:middle,
        ellipsis_line_skip=ellipsis_line_skip,
        equal_columns_width=equal_columns_width,
        title_same_width_as_table=title_same_width_as_table,
        row_label_header_crayon=rlh_crayon,
        title_crayon=title_crayon,
        row_label_crayon=rl_crayon,
        header_crayon=ch_crayon
    )
    if !isnothing(title)
        conf = (conf..., title=title)
    end

    if !isnothing(rows)
        conf = (conf..., row_labels=rows)
    end

    if !isnothing(row_label_title)
        conf = (conf..., row_label_column_title=row_label_title)
    end

    if !isnothing(hlines)
        conf = (conf..., hlines=hlines)
    end



    return conf
end