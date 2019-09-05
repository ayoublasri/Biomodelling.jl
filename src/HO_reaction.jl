function HO_reaction(data)

# Computes the highest order of reaction in which species S_i appears as a
# reatant.


    hor = zeros(data.N, 1)
    mrr = zeros(data.N, 1)

    for j = 1:data.M # run for all reactions

        for i = 1:2 # for the reactant species

            spe_i = data.cm_rea[j, i] # species index

            max_reac_order = maximum(abs.(data.stoichio[j , 1:2])) # maximum order of reactants of reaction j
            reac_order = abs(data.stoichio[j , i]) # order of species i in reaction j

            if spe_i !== 0 && spe_i !== 1 # if i-th species is not zero or NULL

                if hor[spe_i] < max_reac_order # compare the required reactants
                    hor[spe_i] = max_reac_order
                end

                if mrr[spe_i] < reac_order # compare the required reactants
                    mrr[spe_i] = reac_order
                end
            end
        end
    end

    data.hor = hor
    data.mrr = mrr

    return data

end
