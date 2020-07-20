function ssa_switch(data,count)

    temp_data = deepcopy(data)

    temp_species_ts = zeros(temp_data.switch_steps+1, temp_data.N)
    temp_t = zeros(temp_data.switch_steps+1, 1)

    new_species_ts = temp_data.X
    old_species_ts = temp_data.X

    c = zeros(temp_data.M, 1)
    for k = 1:temp_data.M
        c[k] = rates(temp_data, k)
    end
    c_cum = cumsum(c,dims=1)
    c0 = c_cum[end]

    new_t = temp_data.start-temp_data.tau
    old_t = temp_data.start
    i = 0
    for tt = temp_data.start:temp_data.tau:temp_data.switch_steps + temp_data.start - temp_data.tau
        # iterate
        i = i + 1
        while new_t <= tt
            old_t = new_t
            old_species_t = new_species_ts

            if iszero(c0)
                println("System reached a degenerate point. Exiting...")
                species_ts = species_ts[1:i - 1,:]
                t = t[1:i - 1]
            end

            # find which reaction will fire
            reaction = minimum(findall(c0 * rand() .< c_cum))[1] #test

            # compute waiting time
            dt = -log(rand()) / c0

            # update species after reaction
            for k = 1:4
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
            new_species_ts = temp_data.X
            new_t = old_t + dt
        end
        temp_species_ts[i,:] = old_species_ts'
        temp_t[i] = count*temp_data.tau
        count = count + 1
    end
    return temp_t, temp_species_ts
end
