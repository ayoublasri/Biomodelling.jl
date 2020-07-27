function exponential_growth(data,div_noise::Float64,alg::Function)

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
        V1 , V_D, I1, I2= division(V,2.0,output_expression,div_noise,temp_data)
        V, output_expression = replace_cells(V1,V_D,I1,I2)
        output_V[i,:] = V
        output_X[i,:,:] = output_expression
        output_t[i] = (i-1)*temp_data.tau
        expression = output_expression
    end
    return output_t, output_V, output_X
end

function division(V::Array{Float64,1},V_f::Float64,expression,div_noise::Float64,data)
    V_D = []
    expression_D = []
    temp_data = deepcopy(data)
    temp_expression = copy(expression)
    temp_V = copy(V)
    out = V .> V_f
    if any(out)
        cell_div = 0.5.+div_noise.*randn(sum(out.==true))
        V[out] .= cell_div.*temp_V[out]
        exclude = String.(temp_data.species)
        vv = occursin.("on",exclude)
        cc = occursin.("off",exclude)
        aa = vv+cc
        bb=findall(x->x!=1,aa)
        for i in bb
            expression[out,i] .= getBinomial.(expression[out,i],cell_div)
        end
        V_D = zeros(sum(out.==true))
        V_D .= temp_V[out] .- V[out]
        expression_D = zeros(sum(out.==true),size(expression,2))
        expression_D = expression[out,:]
        expression_D = temp_expression[out,bb] .- expression[out,bb]
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
