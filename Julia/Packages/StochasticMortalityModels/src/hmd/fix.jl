export fix_zero_rates!
function fix_zero_rates!(rates::Matrix{Float64})
    issue_idx = findall(r -> r == 0.0, rates)
    for ci in issue_idx
        t = ci.I[2]
        nt = size(rates, 2)
        x = ci.I[1]
        rx = rates[x, :]
        available_rates = findall(x -> x != 0.0, rx)
        forward_idx = (t+1):nt
        backward_idx = (1):(t-1)
        next_idx = findfirst(j -> j in available_rates, forward_idx)
        prev_idx = findlast(j -> j in available_rates, backward_idx)
        try
            if !isnothing(next_idx) && !isnothing(prev_idx)
                next_rate = rates[x, forward_idx[next_idx]]
                prev_rate = rates[x, backward_idx[prev_idx]]
                rates[ci] = 0.5 * (next_rate + prev_rate)
            elseif isnothing(next_idx)
                rates[ci] = rates[x, backward_idx[prev_idx]]
            else
                rates[ci] = rates[x, forward_idx[next_idx]]
            end
        catch error
            println("Got error for x=", x, " t=", t)
            println(error)
        end

    end

    return rates
end