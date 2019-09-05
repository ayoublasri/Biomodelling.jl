function comp_tau(c, data)

    tau = 1e5; # initialize with a big tau

    g = comp_g(data)

    for i = 2:data.N
    # compute h
        hi = max(1, data.epsilon * data.X[i] / g[i])

    # find reactions where species i is involved
        rxns = data.cm_spe[i, 2:1 + data.cm_spe[i,1]]
        rxns1 = unique(rxns)

    # find the respective stoichiometric elements [not optimized]
        idx = findall(x-> x==i, data.cm_rea[rxns1, :])
        tmp_stoichio = data.stoichio[rxns1, :]
        vi = tmp_stoichio[idx]
    # compute mean()
        mi = c[rxns]'* vi
    # compute variance
        sigma2 = c[rxns]' * vi.^2
    # temporary tau

        tau_tmp = min(hi / abs.(mi[1]), hi^2 / sigma2[1])
    # compare minimum(tau) with temporary tau
        if tau_tmp < tau
            tau = tau_tmp
        end
    end

    return tau

end
