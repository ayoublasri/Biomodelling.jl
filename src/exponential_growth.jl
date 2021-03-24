function exponential_growth(data,trans_index::Array{Int64,1},div_noise::Float64,alg::Function,Ni::Float64)

    temp_data = deepcopy(data)

    output_V = zeros(temp_data.NoJ + 1,temp_data.NoC)
    output_X = zeros(Int,temp_data.NoJ +1, temp_data.NoC, temp_data.N)
    output_t = zeros(temp_data.NoJ + 1)
    expression = repeat(temp_data.X',temp_data.NoC)

    exclude = String.(temp_data.species)
    vv = occursin.("on",exclude)
    bb = findall(x->x==1,vv)
    for i in bb
        expression[:,i] .= getBinomial.(ones(Int,temp_data.NoC),0.5)
    end

    cc = occursin.("off",exclude)
    aa = findall(x->x==1,cc)
    expression[:,aa] .= 1 .- expression[:,bb]

    V = ones(temp_data.NoC) .+ (Ni .* ones(temp_data.NoC)) .*rand(temp_data.NoC)
    V_f = 2 .+ 0.001 .* randn(temp_data.NoC)
    output_V[1,:] = V
    output_X[1,:,:] = repeat(temp_data.X',temp_data.NoC)
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        # iterate
        i = i + 1
        output_expression = genexpression(temp_data,expression,trans_index,V,alg::Function)
        V = V .* exp(temp_data.growth_rate*temp_data.tau)
        V1 , V_D, I1, I2= division(V,V_f,output_expression,div_noise,temp_data)
        V, output_expression = replace_cells(V1,V_D,I1,I2)
        output_V[i,:] = V
        output_X[i,:,:] = output_expression
        output_t[i] = (i-1)*temp_data.tau
        expression = output_expression
    end
    return output_t, output_V, output_X
end

function exponential_growth(data,div_noise::Float64,alg::Function,Ni::Float64)

    temp_data = deepcopy(data)
    output_V = zeros(temp_data.NoJ + 1,temp_data.NoC)
    output_X = zeros(Int,temp_data.NoJ +1, temp_data.NoC, temp_data.N)
    output_t = zeros(temp_data.NoJ + 1)
    expression = repeat(temp_data.X',temp_data.NoC)

    exclude = String.(temp_data.species)
    vv = occursin.("on",exclude)
    bb = findall(x->x==1,vv)
    for i in bb
        expression[:,i] .= getBinomial.(ones(Int,temp_data.NoC),0.5)
    end

    cc = occursin.("off",exclude)
    aa = findall(x->x==1,cc)
    expression[:,aa] .= 1 .- expression[:,bb]

    V = ones(temp_data.NoC) .+ (Ni .* ones(temp_data.NoC)) .*rand(temp_data.NoC)
    V_f = 2 .+ 0.001 .* randn(temp_data.NoC)
    output_V[1,:] = V
    output_X[1,:,:] = repeat(temp_data.X',temp_data.NoC)
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        i = i + 1
        output_expression = genexpression(temp_data,expression,alg::Function)
        V = V .* exp(temp_data.growth_rate*temp_data.tau)
        V1 , V_D, I1, I2= division(V,V_f,output_expression,div_noise,temp_data)
        V, output_expression = replace_cells(V1,V_D,I1,I2)
        output_V[i,:] = V
        output_X[i,:,:] = output_expression
        output_t[i] = (i-1)*temp_data.tau
        expression = output_expression
    end
    return output_t, output_V, output_X
end

function exponential_growth_final(data,trans_index::Array{Int64,1},div_noise::Float64,alg::Function,Ni::Float64)

    temp_data = deepcopy(data)

    output_V = zeros(temp_data.NoC)
    output_X = zeros(Int,temp_data.NoC, temp_data.N)
    output_t = zeros(temp_data.NoJ + 1)
    expression = repeat(temp_data.X',temp_data.NoC)

    exclude = String.(temp_data.species)
    vv = occursin.("on",exclude)
    bb = findall(x->x==1,vv)
    for i in bb
        expression[:,i] .= getBinomial.(ones(Int,temp_data.NoC),0.5)
    end

    cc = occursin.("off",exclude)
    aa = findall(x->x==1,cc)
    expression[:,aa] .= 1 .- expression[:,bb]

    V = ones(temp_data.NoC) .+ (Ni .* ones(temp_data.NoC)) .*rand(temp_data.NoC)
    V_f = 2 .+ 0.001 .* randn(temp_data.NoC)
    output_V = V
    output_X = repeat(temp_data.X',temp_data.NoC)
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        # iterate
        i = i + 1
        output_expression = genexpression(temp_data,expression,trans_index,V,alg::Function)
        V = V .* exp(temp_data.growth_rate*temp_data.tau)
        V1 , V_D, I1, I2= division(V,V_f,output_expression,div_noise,temp_data)
        V, output_expression = replace_cells(V1,V_D,I1,I2)
        output_V = V
        output_X = output_expression
        expression = output_expression
    end
    return output_V, output_X
end

function exponential_growth_final(data,div_noise::Float64,alg::Function,Ni::Float64)

    temp_data = deepcopy(data)
    output_V = zeros(temp_data.NoC)
    output_X = zeros(temp_data.NoC, temp_data.N)
    expression = repeat(temp_data.X',temp_data.NoC)

    exclude = String.(temp_data.species)
    vv = occursin.("on",exclude)
    bb = findall(x->x==1,vv)
    for i in bb
        expression[:,i] .= getBinomial.(ones(Int,temp_data.NoC),0.5)
    end

    cc = occursin.("off",exclude)
    aa = findall(x->x==1,cc)
    expression[:,aa] .= 1 .- expression[:,bb]

    V = ones(temp_data.NoC) .+ (Ni .* ones(temp_data.NoC)) .*rand(temp_data.NoC)
    V_f = 2 .+ 0.001 .* randn(temp_data.NoC)
    output_V = V
    output_X = repeat(temp_data.X',temp_data.NoC)
    i = 1
    for tt = temp_data.tau:temp_data.tau:temp_data.T
        i = i + 1
        output_expression = genexpression(temp_data,expression,alg::Function)
        V = V .* exp(temp_data.growth_rate*temp_data.tau)
        V1 , V_D, I1, I2= division(V,V_f,output_expression,div_noise,temp_data)
        V, output_expression = replace_cells(V1,V_D,I1,I2)
        output_V = V
        output_X = output_expression
        expression = output_expression
    end
    return output_V, output_X
end

function division(V::Array{Float64,1},V_f::Array{Float64,1},expression,div_noise::Float64,temp_data)
    V_D = []
    expression_D = []
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
        bb = findall(x->x!=1,aa)
        for i in bb
            expression[out,i] .= getBinomial.(expression[out,i],cell_div)
        end
        V_D = zeros(sum(out.==true))
        V_D .= temp_V[out] .- V[out]
        expression_D = zeros(sum(out.==true),size(expression,2))
        expression_D = expression[out,:]
        expression_D[:,bb] = temp_expression[out,bb] .- expression[out,bb]
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

function genexpression(temp_data,expression,trans_index::Array{Int64,1},V::Array{Float64,1},alg::Function)
    temp_data.T = temp_data.tau
    temp_data.NoJ = 1
    output_expression = zeros(Int,temp_data.NoC,temp_data.N)
    trans_rates = temp_data.kr[trans_index,1]
    for i = 1:temp_data.NoC
        temp_data.kr[trans_index,1] = V[i] * trans_rates
        temp_data.X[1:end] = expression[i,:]
        if alg == Biomodelling.tauleapswitch || alg == Biomodelling.non_negative_Poisson_tauleap
            temps, valeur = alg(temp_data,1)
        elseif alg == Biomodelling.tauleapswitch2
            valeur = alg(temp_data,1)
        elseif alg == Biomodelling.ssa3 
            valeur = alg(temp_data)
        else
            temps, valeur = alg(temp_data)
        end
        for j = 1:temp_data.N
            if alg == Biomodelling.ssa3 || alg == Biomodelling.tauleapswitch2
                output_expression[i,j] = valeur[j]
            else
                output_expression[i,j] = valeur[end,j]
            end
        end
    end
    temp_data.kr[trans_index,1] = trans_rates
    return output_expression
end

function genexpression(temp_data,expression,alg::Function)
    temp_data.T = temp_data.tau
    temp_data.NoJ = 1
    output_expression = zeros(Int,temp_data.NoC,temp_data.N)
    for i = 1:temp_data.NoC
        temp_data.X[1:end] = expression[i,:]
        if alg == Biomodelling.tauleapswitch ||  alg == Biomodelling.non_negative_Poisson_tauleap
            temps, valeur = alg(temp_data,1)
        elseif alg == Biomodelling.tauleapswitch2
            valeur = alg(temp_data,1)
        elseif alg == Biomodelling.ssa3 
            valeur = alg(temp_data)
        else
            temps, valeur = alg(temp_data)
        end
        for j = 1:temp_data.N
            if alg == Biomodelling.ssa3 || alg == tauleapswitch2
                output_expression[i,j] = valeur[j]
            else
                output_expression[i,j] = valeur[end,j]
            end
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
