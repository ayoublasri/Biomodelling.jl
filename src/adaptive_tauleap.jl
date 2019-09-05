function adaptive_tauleap(data)

    temp_data = deepcopy(data)

    # ---- OUTPUT PARAMETERS ----
    #   output.t: time vector, [NoJx1]
    #   output.species_ts: the concentrations of each species at each time step, [NoJxN]


    output_species_ts = zeros(temp_data.NoJ+1, temp_data.N)
    output_t = zeros(temp_data.NoJ+1, 1)
    species_ts = zeros(size(temp_data.X))
    output_species_ts[1,:] = temp_data.X
    output_t[1] = 0

    # Rates initialization
    c = zeros(temp_data.M, 1)
    for k = 1:temp_data.M
        c[k] = rates(temp_data, k)
    end
    t = 0.0
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        # iterate
        i = i + 1
        while t < tt
            # check if rates are NaN()
            if any(isnan, c)
                error("Rates in simulation became Nan.")
            end

            # check if there is nowhere to move()
            if iszero(sum(c))
                println("System reached a degenerate point. Exiting...");
                species_ts = species_ts[1:i - 1,:]
                t = t[1:i - 1]
                return t
            end

            # compute tau from exponential
            dt = comp_tau(c, temp_data)
            if t + dt > tt
                dt = tt-t
            end

            # find which reaction will fire
            r = [PoissonRandom.pois_rand(v * temp_data.tau) for v in c]

            # update species after reaction
            for reaction = 1:temp_data.M
                if !iszero(r[reaction])
                    for k = 1:4
                        ind = temp_data.cm_rea[reaction,k]
                        if !iszero(temp_data.cm_rea[reaction,k])
                            temp_data.X[ind] = temp_data.X[ind] + r[reaction] * temp_data.stoichio[reaction,k]
                        end
                    end
                end
            end


            # update rates if there was a reaction fired
            for reaction = 1:temp_data.M
                if !iszero(r[reaction])
                    c = update_rates(c, temp_data, reaction)
                end
            end
            # update output
            species_ts = temp_data.X
            t = t + dt
            if t == tt
                output_species_ts[i,:] = species_ts
                output_t[i] = t
            end
        end
    end
    return output_t, output_species_ts
end
