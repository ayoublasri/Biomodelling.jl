function ssa_switch2(data,count)

    temp_data = deepcopy(data)

    species_ts = zeros(size(temp_data.X))

    c = zeros(temp_data.M, 1)
    for k = 1:temp_data.M
        c[k] = rates(temp_data, k)
    end
    c_cum = cumsum(c,dims=1)
    c0 = c_cum[end]

    t = temp_data.start
    i = 0
    while i <= count
        i = i + 1
        if iszero(c0)
            println("System reached a degenerate point. Exiting...")
            species_ts = species_ts[1:i - 1,:]
            t = t[1:i - 1]
        end

        # find which reaction will fire
        reaction = minimum(findall(c0 * rand() .< c_cum))[1]

        # compute waiting time
        dt = -log(rand()) / c0

        # update species after reaction
        for k = 1:10
            ind = temp_data.cm_rea[reaction, k]
            if !iszero(temp_data.cm_rea[reaction, k])
                temp_data.X[ind] = temp_data.X[ind] + temp_data.stoichio[reaction, k]
            end
        end

        # update rates involving species of 'reaction' after species update
        c = update_rates(c, temp_data, reaction)
        c_cum = cumsum(c,dims=1)
        c0 = c_cum[end]

        # update output
        species_ts = temp_data.X
        t = t + dt
    end
    return t, species_ts
end
