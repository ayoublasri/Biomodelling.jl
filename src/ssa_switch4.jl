function ssa_switch4(temp_data,count)

    new_species_ts = temp_data.X
    old_species_ts = temp_data.X

    temp_t = zeros(temp_data.switch_steps, 1)

    c = zeros(temp_data.M, 1)
    for k = 1:temp_data.M
        c[k] = rates(temp_data, k)
    end
    c_cum = cumsum(c,dims=1)
    c0 = c_cum[end]

    new_t = temp_data.start
    old_t = temp_data.start
    i = 0
    for tt = temp_data.start:temp_data.tau:temp_data.switch_steps*temp_data.tau + temp_data.start-temp_data.tau
        # iterate
        i = i + 1
        while new_t <= tt
            old_t = new_t
            old_species_ts = new_species_ts

            if iszero(c0)
                println("System reached a degenerate point. Exiting...")
                species_ts = species_ts[1:i - 1,:]
                t = t[1:i - 1]
            end

            # find which reaction will fire
            drxn = rand() * c0
            psm = 0.0
            reaction = 0
            while psm < drxn
                reaction = reaction + 1
                psm = psm + c[reaction]
            end
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

        temp_t = count*temp_data.tau
        count = count + 1
    end
    return temp_t, old_species_ts
end
