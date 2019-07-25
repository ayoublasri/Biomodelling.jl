function rates(data, reaction)

    rate = data.kr[reaction]

    #if reaction == 2 # 2nd reaction has a special nonlinear rate. It should
    # hold the correspodence: p53==X(2) & mdm2==X(4)
    #rate = data.kr[reaction] + 1.7 * data.X[4] / (data.X[2] + 0.01)
    #else
    #rate = data.kr[reaction]
    #end

    for j=1:2
        species = data.cm_rea[reaction, j]

        # check if species is not NULL [species==1]
        if species > 1

            # check if there are enough elements of species
            if data.X[species] < abs(data.stoichio2[reaction, j])
                rate = 0
            else
                if abs(data.stoichio2[reaction,j])==0
                    rate = rate
                elseif abs(data.stoichio2[reaction,j])==1
                    rate = rate * data.X[species]
                elseif abs(data.stoichio2[reaction,j])==2
                    rate = rate * data.X[species] * (data.X[species]-1) / 2
                elseif abs(data.stoichio2[reaction,j])==3
                    rate = rate * data.X[species] *(data.X[species]-1) * (data.X[species]-2) / 6
                else
                    rate = rate*binomial(data.X[species],abs(data.stoichio2[reaction,j]))
                end
            end
        end

    end

    return rate

end
