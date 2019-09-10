function exponential_growth(data)

    temp_data = deepcopy(data)

    output_V = zeros(temp_data.NoJ + 1,temp_data.NoC)
    output_t = zeros(temp_data.NoJ + 1)

    V = ones(temp_data.NoC)
    output_V[1,:] = V
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        # iterate
        i = i + 1
        V = V .*exp(temp_data.growth_rate*temp_data.tau)
        V1 , V_D = division(V,2.0)
        V = replace_cells(V1,V_D)
        output_V[i,:] = V
        output_t[i] = (i-1)*temp_data.tau
    end
    return output_t, output_V
end

function division(V::Array{Float64,1},V_f::Float64)
    V_D = []
    temp_V = copy(V)
    out = V .> V_f
    if any(out)
        V[out].=0.5.*temp_V[out]
        V_D = zeros(length(out))
        V_D .= temp_V[out] .- V[out]
    end
    return V, V_D
end

function replace_cells(V,V_D)
    if !isempty(V_D)
        V[Random.randperm(length(V))[1:length(V_D)]] .= V_D
    end
    return V
end
