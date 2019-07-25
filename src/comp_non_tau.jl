function comp_non_tau(c, data, critical_reaction)

                g = comp_g(data)
                tau = data.tau

                non_critical_reaction = Array{Int64,1}()
                temp_non_critical = collect(1:data.M)
                for i = 1:length(critical_reaction)
                    temp_non_critical = filter(!x->x==critical_reaction[i],temp_non_critical)
                end
                non_critical_reaction = temp_non_critical

                for i = 2:data.N
                                # compute h
                                hi = max(1, data.epsilon * data.X[i] / g[i])

                                # find non critical reactions where species i is involved
                                rxns = non_critical_reaction
                                # find the respective stoichiometric elements [not optimized]
                                idx = findall(x-> x==i, data.cm_rea[rxns, :])
                                if length(idx)!==0
                                tmp_stoichio = data.stoichio[rxns, :]
                                vi = tmp_stoichio[idx]
                                # compute mean()
                                mi = c[rxns]'* vi
                                # compute variance
                                sigma2 = c[rxns]'* vi.^2
                                # temporary tau
                                tau_tmp = min(hi / abs.(mi[1]), hi^2 / sigma2[1])
                                # compare minimum(tau) with temporary tau
                                if tau_tmp < data.tau
                                                tau = tau_tmp
                                end
                                end
                end
                return tau
end
