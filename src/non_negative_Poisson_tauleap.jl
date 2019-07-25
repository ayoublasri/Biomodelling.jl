function non_negative_Poisson_tauleap(data,n_c)

    temp_data = deepcopy(data)

    output_species_ts = zeros(temp_data.NoJ+1, temp_data.N)
    output_t = zeros(temp_data.NoJ+1, 1)

    species_ts = temp_data.X
    output_species_ts[1,:] = temp_data.X
    output_t[1] = 0

    c = zeros(temp_data.M, 1)
    for k = 1:temp_data.M
        c[k] = rates(temp_data, k)
    end
    # compute critical reactions
    t = 0.0
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        # iterate
        i = i + 1
        while t < tt
            L = compute_L(temp_data)

            indx = findall(x->x<n_c,L)

            critical_reaction = Array{Int64,1}()

            w = Array{Int64,1}()

            for j = 1:length(indx)
                if c[indx[j][1],1] > 0
                    append!(critical_reaction,indx[j][1])
                end
            end
            if length(critical_reaction)== 0
                dt = comp_tau(c, temp_data)
                r = [rand(Distributions.Poisson(v * dt)) for v in c]

                if t + dt > tt
                    dt = tt-t
                end

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
                for reaction = 1:temp_data.M
                    if !iszero(r[reaction])
                        c = update_rates(c, temp_data, reaction)
                    end
                end
                species_ts = temp_data.X
                t = t + dt
            else
                dt = comp_non_tau(c,temp_data,critical_reaction)
                if dt < 10*(1/sum(c))

                    temp_data.start = t
                    t_out, X_out = ssa_switch2(temp_data,10)
                    temp_data.X = X_out
                    t = t_out
                    species_ts = temp_data.X

                else
                    for i in critical_reaction
                        append!(w,c[i])
                    end
                    c_critical = cumsum(w,dims=1)
                    dt2 = -log(rand()) / c_critical[end]

                    if dt < dt2
                        temp_data.tau = dt
                        r = [rand(Distributions.Poisson(v * temp_data.tau)) for v in c]
                        for i in critical_reaction
                            r[i] = 0
                        end
                    elseif dt2 <= dt
                        temp_data.tau = dt2
                        r = [rand(Distributions.Poisson(v * temp_data.tau)) for v in c]
                        reaction = minimum(findall(c_critical[end] * rand() .< c_critical))[1]
                        for i in critical_reaction
                            if i == reaction
                                r[i] = 1
                            else
                                r[i] = 0
                            end
                        end
                    end

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
                    for reaction = 1:temp_data.M
                        if !iszero(r[reaction])
                            c = update_rates(c, temp_data, reaction)
                        end
                    end
                    species_ts = temp_data.X
                    t = t + dt
                end
            end
        end
        # update output
        if t >= tt
            output_species_ts[i,:] = species_ts
            output_t[i] = tt
        end
    end
    return output_t, output_species_ts
end
