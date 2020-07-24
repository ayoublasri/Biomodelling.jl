function ssa(data)

    temp_data = deepcopy(data)

    output_species_ts = zeros(Int,temp_data.NoJ + 1, temp_data.N)
    output_t = zeros(temp_data.NoJ + 1)

    new_species_ts = temp_data.X
    old_species_ts = temp_data.X
    output_species_ts[1,:] = temp_data.X
    output_t[1] = 0

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
            old_species_t = new_species_ts
            if iszero(c0)
                println("System reached a degenerate point. Exiting...")
                output_species_ts = output_species_ts[1:i - 1,:]
                output_t = output_t[1:i - 1]
            end
            
            # find which reaction will fire
            reaction = minimum(findall(c0 * rand() .< c_cum))[1]
            # compute waiting time
            dt = -log(rand())/c0
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
        output_species_ts[i,:] = old_species_ts
        output_t[i] = (i-1)*temp_data.tau
    end
    return output_t, output_species_ts
end
