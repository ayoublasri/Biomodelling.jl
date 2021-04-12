function rates(data, reaction)

    rate = data.kr[reaction,:]

    if contains(data.model[reaction].name, "comb_a")
        reac = data.cm_rea[reaction,1:length(data.model[reaction].reactants)]
        reac = reac .- size(data.matrix,1)
        produit = zeros(length(reac))
        for i = 1:length(reac)
            produit[i] = ((data.X[reac[i]]/rate[4])^(rate[2])) / (rate[3]^(rate[2]) + (data.X[reac[i]]/rate[4])^(rate[2]))
        end
        rate = rate[1] * rate[4] * prod(produit)
        return rate
    elseif contains(data.model[reaction].name, "comb_i")
        reac = data.cm_rea[reaction,1:length(data.model[reaction].reactants)]
        reac = reac .- size(data.matrix,1)
        produit = zeros(length(reac))
        for i = 1:length(reac)
            produit[i] = (rate[3]^(rate[2])) / (rate[3]^(rate[2]) + (data.X[reac[i]]/rate[4])^(rate[2]))
        end
        rate = rate[1] * rate[4] * prod(produit)
        return rate
    elseif contains(data.model[reaction].name, "comb_rea")
        reac = data.cm_rea[reaction,1:length(data.model[reaction].reactants)]
        reac = reac .- size(data.matrix,1)
        dim1 = data.cm_rea[reaction,6] - 1
        dim2 = reac .- 1
        produit = zeros(length(dim2))
        for i = 1:length(dim2)
            if data.matrix[dim2[i],dim1] == 1
                produit[i] = ((data.X[reac[i]]/rate[4])^(rate[2])) / (rate[3]^(rate[2]) + (data.X[reac[i]]/rate[4])^(rate[2]))
            elseif  data.matrix[dim2[i],dim1] == -1
                produit[i] = (rate[3]^(rate[2])) / (rate[3]^(rate[2]) + (data.X[reac[i]]/rate[4])^(rate[2]))
            end
        end
        rate = rate[1] * rate[4] * prod(produit)
        return rate
    else
        for j = 1:5
            species = data.cm_rea[reaction, j]

            # check if species is not NULL [species==1]
            if species > 0

                # check if there are enough elements of species
                if data.X[species] < abs(data.stoichio2[reaction, j])
                    rate = 0
                else
                    if contains(data.model[reaction].name, "act") 
                        rate = rate[1] * rate[4] * ((data.X[species]/rate[4])^(rate[2])) / (rate[3]^(rate[2]) + (data.X[species]/rate[4])^(rate[2]))
                    elseif contains(data.model[reaction].name, "inhib") 
                        rate = rate[1] * rate[4] * (rate[3]^(rate[2])) / (rate[3]^(rate[2]) + (data.X[species]/rate[4])^(rate[2]))
                    else
                        if abs(data.stoichio2[reaction,j]) == 0
                            rate = rate[1]
                        elseif abs(data.stoichio2[reaction,j]) == 1
                            rate = rate[1] * data.X[species]
                        elseif abs(data.stoichio2[reaction,j]) == 2
                            rate = rate[1] * data.X[species] * (data.X[species] - 1) / 2
                        elseif abs(data.stoichio2[reaction,j]) == 3
                            rate = rate[1] * data.X[species] * (data.X[species] - 1) * (data.X[species] - 2) / 6
                        else
                            rate = rate[1] * binomial(data.X[species], abs(data.stoichio2[reaction,j]))
                        end
                    end
                end
                return rate
            end
        end
    end
end
