function tauleapswitch2(temp_data,ssa_steps)
    #temp_data = deepcopy(data)
    temp_X = temp_data.X

    # Rates initialization
    c = zeros(temp_data.M, 1)
    for reaction = 1:temp_data.M
        c[reaction] = rates(temp_data, reaction)
    end
    i = 2
    t = 0.0
    species_ts = temp_data.X
    # iterate
    while i <= temp_data.NoJ +1
        # check if rates are NaN()
        if any(isnan, c)
            throw(error("Rates in simulation became Nan."))
        end

        # check if there is nowhere to move()
        if iszero(sum(c))
            println("System reached a degenerate point. Exiting...");
        end
        # find which reaction will fire ()
        r = [PoissonRandom.pois_rand(v * temp_data.tau) for v in c]

        # update species after reaction
        for reaction = 1:temp_data.M
            if !iszero(r[reaction])                # if is not zero we update the reaction
                for k = 1:4
                    ind = temp_data.cm_rea[reaction,k]
                    if !iszero(temp_data.cm_rea[reaction,k])
                        temp_X[ind] = temp_data.X[ind] + r[reaction] * temp_data.stoichio[reaction,k]
                    end
                end
            end
        end

        if any(temp_X .< 0)
            count = i - 1
            temp_data.start = t
            temp_data.switch_steps = min(ssa_steps,temp_data.NoJ-count+1)
            temp_data.X[1:end] = species_ts

            output_t, output_species_ts=ssa_switch4(temp_data,count)
            t = output_t
            species_ts = output_species_ts
            temp_data.X[1:end] = species_ts
            i = i + temp_data.switch_steps
            if t > temp_data.T
                return species_ts
            end
        else
            temp_data.X = temp_X
            for reaction = 1:temp_data.M
                if !iszero(r[reaction]) # if there was a reaction fired
                    c = update_rates(c, temp_data, reaction)
                end
            end

            # update output
            species_ts = temp_data.X
            t =  t + temp_data.tau

            if t > temp_data.T
                return species_ts
            end
            i = i+1
        end
    end
    return species_ts
end
