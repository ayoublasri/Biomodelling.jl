function exponential_growth(data,alg::Function)

    temp_data = deepcopy(data)

    output_V = zeros(temp_data.NoJ + 1,temp_data.NoC)
    output_X = zeros(Int,temp_data.NoJ +1, temp_data.NoC, temp_data.N)
    output_t = zeros(temp_data.NoJ + 1)
    expression = repeat(temp_data.X',temp_data.NoC)

    V = ones(temp_data.NoC)
    output_V[1,:] = V
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        # iterate
        i = i + 1
        output_expression = genexpression(temp_data,expression,alg::Function)
        V = V .*exp(temp_data.growth_rate*temp_data.tau)
        V1 , V_D, I1, I2= division(V,2.0,output_expression)
        V, output_expression = replace_cells(V1,V_D,I1,I2)
        output_V[i,:] = V
        output_X[i,:,:] = output_expression
        output_t[i] = (i-1)*temp_data.tau
        expression = output_expression
    end
    return output_t, output_V, output_X
end

function division(V::Array{Float64,1},V_f::Float64,expression)
    V_D = []
    expression_D = []
    temp_expression = copy(expression)
    temp_V = copy(V)
    out = V .> V_f
    if any(out)
        V[out].=0.5.*temp_V[out]
        for i = 1: size(expression,2)
            expression[out,i] .= getBinomial.(expression[out,i],0.5)
        end
        V_D = zeros(length(out))
        V_D .= temp_V[out] .- V[out]
        expression_D = zeros(size(expression[out,:]))
        expression_D = temp_expression[out,:] .- expression[out,:]
    end
    return V, V_D, expression, expression_D
end

function replace_cells(V,V_D,expression,expression_D)
    if !isempty(V_D)
        indexes = Random.randperm(length(V))[1:length(V_D)]
        V[indexes] .= V_D
        expression[indexes,:] .= expression_D
    end
    return V, expression
end

function genexpression(data,expression,alg::Function)
    temp_data = deepcopy(data)
    temp_data.T = temp_data.tau
    temp_data.NoJ = 1
    output_expression = zeros(Int,temp_data.NoC,temp_data.N)
    for i = 1:temp_data.NoC
        temp_data.X[:,1] = expression[i,:]
        if alg == Biomodelling.tauleapswitch
            temps, valeur = alg(temp_data,100)
        elseif alg == Biomodelling.non_negative_Poisson_tauleap
            temps, valeur = alg(temp_data,10)
        else
            temps, valeur = alg(temp_data)
        end
        for j = 1:temp_data.N
            output_expression[i,j] = valeur[end,j]
        end
    end
    return output_expression
end

function getBinomial(n::Integer, p::Real)
   log_q = log(1.0 - p)
   x = 0
   sum = 0.0
   while true
       sum += log(rand()) / (n - x)
       sum < log_q && break
       x += 1
   end
   return x
end
