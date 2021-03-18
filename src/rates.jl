function rates(data, reaction)

    rate = data.kr[reaction,:]

    for j = 1:2
        species = data.cm_rea[reaction, j]

        # check if species is not NULL [species==1]
        if species > 0

            # check if there are enough elements of species
            if data.X[species] < abs(data.stoichio2[reaction, j])
                rate = 0
            else
                if contains(data.model[reaction].name, "act") 
                    rate = rate[1] * data.X[species]^(rate[2]) / (rate[3]^(rate[2]) + data.X[species]^(rate[2]))
                elseif contains(data.model[reaction].name, "inhib") 
                    rate = rate[1] * rate[3]^(rate[2]) / (rate[3]^(rate[2]) + data.X[species]^(rate[2]))
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

        end
        return rate

    end
end
