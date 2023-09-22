export readfile

function readfile(
    path::AbstractString;
    remove_missing::Bool=true
)::DataFrame

    if startswith(path, "https")
        println("Downloading file at url `", path, "`")
        data, headers = DelimitedFiles.readdlm(download(path); header=true, skipstart=2)
    else
        data, headers = DelimitedFiles.readdlm(path; header=true, skipstart=2)
    end

    data_dict::Dict{String,Any} = Dict()

    data[:, 2] = ([typeof(age) <: Number ? age : occursin("-", age) ? parse(Float64, split(age, "-")[1]) : age for age in data[:, 2]])
    data[:, 2] = ([age == "110+" ? 110 : age for age in data[:, 2]])

    data[:, 1] = (data[:, 1])

    nrow = size(data, 1)
    select_idx = repeat([true], nrow)
    if remove_missing
        for i in 1:nrow
            if any(==("."), data[i, 3:end])
                select_idx[i] = false
            end
        end
    else
        for i in 3:lastindex(headers[1, :])
            data[:, i] = map(x -> (x == ".") ? missing : Float64(x), data[:, i])
        end
    end

    for c in 1:lastindex(headers[1, :])
        col = data[select_idx, c]
        heading = headers[1, c]
        if (c == 1 || c == 2)
            data_dict[heading] = Vector{Int}(col)
        elseif remove_missing
            data_dict[heading] = Vector{Real}(col)
        else
            data_dict[heading] = Vector{Union{Real,Missing}}(col)
        end
    end


    df = DataFrame(data_dict)

    col_order = Vector{Int64}(indexin(headers, names(df))[1, :])

    df = df[:, col_order]

    return df
end