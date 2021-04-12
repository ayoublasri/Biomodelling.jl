function compute_L(data)
    temp_L = Array{Float64,1}()
    L = zeros(data.M,1)
    for i =1:data.M
        if any(data.stoichio2[i,1:5].<0.0)
            ind = findall(data.stoichio2[i,1:5].<0.0)
            for k = 1:length(ind)
                append!(temp_L,round(data.X[data.cm_rea[i,ind[k]]]/abs(data.stoichio2[i,ind[k]])))
            end
            L[i] = minimum(temp_L)
        else
            L[i] = 1000
        end
    end
    return L
end
