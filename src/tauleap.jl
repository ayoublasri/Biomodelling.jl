function tauleap(data)
    temp_data = deepcopy(data)
# Stochastic tau-leap algorithm for well-mixed reaction systems.

# ---- OUTPUT PARAMETERS ----
#   output.t: time vector, [NoJx1]
#   output.species_ts: the concentrations of each species at each time step, [NoJxN]


# initialize output struct()
    species_ts = zeros(temp_data.NoJ+1, temp_data.N)
    t          = zeros(temp_data.NoJ+1, 1)

# initial species population and initial time
    species_ts[1,:] =  temp_data.X
    t[1]            =  0

# Rates initialization
    c = zeros(temp_data.M, 1)
    for reaction = 1:temp_data.M
        c[reaction] = rates(temp_data, reaction)
    end

# iterate
    for i = 2:temp_data.NoJ+1
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
                        temp_data.X[ind] = temp_data.X[ind] + r[reaction] * temp_data.stoichio[reaction,k]
                        if temp_data.X[ind] < 0
                            species_ts = species_ts[1:i,:]
                            t = t[1:i]
                            error("Ended at $i. Try smaller tau")
                        end
                    end
                end
            end
        end
    # update rates

        for reaction = 1:temp_data.M
            if !iszero(r[reaction]) # if there was a reaction fired
                c = update_rates(c, temp_data, reaction)
            end
        end

    # update output

        species_ts[i,:]  =  temp_data.X
        t[i]             =  t[i - 1] + temp_data.tau

        if t[i] > temp_data.T
            species_ts = species_ts[1:i,:]
            t = t[1:i]
        end
    end
    return t, species_ts
end