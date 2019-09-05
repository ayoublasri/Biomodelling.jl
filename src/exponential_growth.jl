function exponential_growth(data)

    temp_data = deepcopy(data)

    output_V = zeros(temp_data.NoJ + 1,temp_data.NoC)
    output_t = zeros(temp_data.NoJ + 1)

    V = ones(temp_data.NoC)
    new_t = 0.0
    old_t = 0.0
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        # iterate
        i = i + 1
        for j in temp_data.NoC
            V[j] = V[j]*exp(temp_data.growth*tt)
            if V[j] > 2.0
                temp_V = 0.5*V[j]
                V[j] = temp_V
                V[rand(1:temp_data.NoC)] = temp_V
            end
        end
        output_V[i,:] = V
        output_t[i] = (i-1)*temp_data.tau
    end
    return output_t, output_V
end
