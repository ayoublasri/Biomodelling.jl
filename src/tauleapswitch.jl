function tauleapswitch(data,ssa_steps)
    temp_data = deepcopy(data)
    temp_X = temp_data.X
    # Stochastic tau-leap algorithm for well-mixed reaction systems.

    # ---- OUTPUT PARAMETERS ----
    #   output.t: time vector, [NoJx1]
    #   output.species_ts: the concentrations of each species at each time step, [NoJxN]


    # initialize output struct()
    species_ts = zeros(temp_data.NoJ +1, temp_data.N)
    t          = zeros(temp_data.NoJ +1, 1)

    # initial species population and initial time
    species_ts[1,:] =  temp_data.X
    t[1]            =  0


    # Rates initialization
    c = zeros(temp_data.M, 1)
    for reaction = 1:temp_data.M
        c[reaction] = rates(temp_data, reaction)
    end
    i = 2
    # iterate
    while i <= temp_data.NoJ +1
        # check if rates are NaN()
        if any(isnan, c)
            throw(error("Rates in simulation became Nan."))
        end

        # check if there is nowhere to move()
        if iszero(sum(c))
            println("System reached a degenerate point. Exiting...");
            species_ts = species_ts[1:i - 1,:]
            t = t[1:i - 1]
            return t
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
            t_tau = t[1:count]
            temp_data.start = t_tau[end]
            temp_data.switch_steps = min(ssa_steps,temp_data.NoJ-count+1)

            temp_data.X = species_ts[count,:]'
            species_tau = species_ts[1:count,1:temp_data.N]

            output_t, output_species_ts=ssa_switch(temp_data,count)
            t[count+1:temp_data.switch_steps+count] = output_t'
            species_ts[count+1: temp_data.switch_steps + count, 1:temp_data.N] = output_species_ts
            temp_data.X = species_ts[end,1:temp_data.N]'
            i = i + temp_data.switch_steps
            if i > temp_data.NoJ +1
                return t, species_ts
            end
        else
            temp_data.X = temp_X
            for reaction = 1:temp_data.M
                if !iszero(r[reaction]) # if there was a reaction fired
                    c = update_rates(c, temp_data, reaction)
                end
            end

            # update output
            species_ts[i,:]  =  temp_data.X
            t[i]             =  t[i-1] + temp_data.tau

            if t[i] > temp_data.T
                species_ts = species_ts[1:i,:]
                t = collect((0:1:i).*temp_data.tau)
            end
            i = i+1
        end
    end
    return t, species_ts
end
