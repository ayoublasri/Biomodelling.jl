function ssa_switch3(data)

    temp_data = deepcopy(data)

    species_ts = zeros(size(temp_data.X))
    output_species_ts = zeros(temp_data.NoJ + 1, temp_data.N)
    output_t = zeros(temp_data.NoJ + 1, 1)

    c = zeros(temp_data.M, 1)
    for k = 1:temp_data.M
        c[k] = rates(temp_data, k)
    end
    c_cum = cumsum(c,dims=1)
    c0 = c_cum[end]

    t = temp_data.start
    i = 0
    while t <= data.T
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
        species_ts = temp_data.X
        t = t + dt
        if i > temp_data.NoJ
            output_species_ts = [output_species_ts;species_ts']
            output_t = [output_t;t]
        end
        if t >= data.T
            itp1 = interpolate(output_species_ts[:,2], BSpline(Linear()))
            itp2 = interpolate(output_t, BSpline(Linear()))
            protein = itp1(range(0,stop=temp_data.T,length=temp_data.NoJ))
            time = itp2(range(0,stop=temp_data.T,length=temp_data.NoJ))
            return time, protein
        end
    end
end
