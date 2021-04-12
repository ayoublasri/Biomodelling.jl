function ssa3(temp_data)

    new_species_ts = temp_data.X

    c = zeros(temp_data.M, 1)
    for k = 1:temp_data.M
        c[k] = rates(temp_data, k)
    end
    c_cum = cumsum(c,dims=1)
    c0 = c_cum[end]
    new_t = 0.0
    old_t = 0.0
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        # iterate
        i = i + 1
        while new_t < tt
            old_t = new_t
            if iszero(c0)
                println("System reached a degenerate point. Exiting...")
            end

            # find which reaction will fire
            # reaction = minimum(findall(c0 * rand() .< c_cum))[1]
            # compute waiting time
            drxn = rand() * c0
            psm = 0.0
            reaction = 0
            while psm < drxn
                reaction = reaction + 1
                psm = psm + c[reaction]
            end
            dt = -log(rand())/c0
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
            new_species_ts = temp_data.X[1:end]
            new_t = old_t + dt
        end
    end
    return new_species_ts
end