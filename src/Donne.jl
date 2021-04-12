mutable struct Donne

    M::Int64
    N::Int64
    T::Float64
    tau::Float64
    start::Float64
    switch_steps::Int64
    epsilon::Float64
    NoJ::Int64
    NoC::Int64
    growth_rate::Float64
    species::Array{Symbol,1}
    X::Matrix{Int64}
    stoichio::Matrix{Int64}
    stoichio2::Matrix{Int64}
    kr::Matrix{Float64}
    cm_rea::Matrix{Int64}
    cm_spe::Matrix{Int64}
    hor::Matrix{Int64}
    mrr::Matrix{Int64}
    model:: Any
    matrix::Matrix{Int64}

    #   data.M: number of reactions, [1x1]
    #   data.N: number of species, [1x1]
    #   data.X: initial population of the system(), [Nx1]
    #   data.stoichio: stoichiometric matrix, [Mx4]
    #   data.kr: rate constants, [Mx1]
    #   data.delays: reaction delays [type of delay and uniform dirstibution parameters],[Mx3]
    #   data.cm_rea: connectivity matrix for the reactions, [Mx4]
    #   data.cm_spe: connectivity matrix for the species, [MxK+1]

    function Donne(model,matrix, initiale, T, tau, NoC, growth_data, epsilon=0.03)

        M = length(model)
        species = getSpecies(model)
        #####################################
        #### load stoichoimetric matrix #####
        #####################################

        N = length(species)

        tmp_stoichio = getStoichio(model)
        tmp_stoichio2 = getStoichio2(model)

        stoichio = zeros(M, 10)
        stoichio2 = zeros(M,10)

        for m = 1:M
            stoichio[tmp_stoichio[m,1], :] = tmp_stoichio[m, 2:11]
        end
        for m = 1:M
            stoichio2[tmp_stoichio2[m,1], :] = tmp_stoichio2[m, 2:11]
        end
        #######################################
        #### load reaction rate constants #####
        #######################################

        kr = zeros(M, 4)
        for i = 1:M
            kr[i,1:3] .= model[i].rate
        end

        ############################################
        #### load reaction connectivity matrix #####
        ############################################
        species = getSpecies(model)
        tmp_cm_rea = getConnectivity(model,species)

        cm_rea = zeros(Int64, M, 10)
        for m = 1:M
            cm_rea[tmp_cm_rea[m,1], :] = tmp_cm_rea[m, 2:11]
        end


        ##################################
        #### load initial population #####
        ##################################
        tmp_X = getInitiale(initiale,species)

        X = zeros(N, 1)
        for j = 1:N
            if tmp_X[j,1] > N
                println("Initial population input file is wrong!")
            end
            X[tmp_X[j,1]] = tmp_X[j, 2]
        end

        #############################################
        #### create species connectivity matrix #####
        #############################################
        K = 50

        cm_spe = zeros(Int64, N, K)
        no_rea = zeros(Int64, N, 1)
        for m = 1:M
            for k = 1:10
                if cm_rea[m,k] !== 0 && cm_rea[m,k] !== 1

                    cm_spe[cm_rea[m,k], no_rea[cm_rea[m,k]] + 2] = m

                    no_rea[cm_rea[m,k]] = no_rea[cm_rea[m,k]] + 1

                    if no_rea[cm_rea[m,k]] > K
                        cm_spe = [cm_spe, zeros(N, K)]
                    end
                end

            end
        end

        cm_spe[:,1] = no_rea

        maxK = maximum(no_rea)
        cm_spe = cm_spe[:,1:maxK + 1]

        NoJ =  Int(T / tau)
        start = tau
        switch_steps = 10
        hor = zeros(Int64, N, 1)
        mrr = zeros(Int64, N, 1)

        for j = 1:M 

            for i = 1:5 

                spe_i = cm_rea[j, i] 

                max_reac_order = maximum(abs.(stoichio[j , 1:5])) 
                reac_order = abs(stoichio[j , i]) 

                if spe_i !== 0 && spe_i !== 1 

                    if hor[spe_i] < max_reac_order
                        hor[spe_i] = max_reac_order
                    end

                    if mrr[spe_i] < reac_order 
                        mrr[spe_i] = reac_order
                    end
                end
            end
        end
        if typeof(growth_data) == Float64
            growth_rate = growth_data
        else
            growth_rate = growth_estimate(growth_data)
        end
        new(M, N, T, tau, start, switch_steps, epsilon, NoJ, NoC, growth_rate, species, X, stoichio, stoichio2,kr, cm_rea, cm_spe, hor, mrr,model,matrix)
    end
end


function getSpecies(model)
    n = length(model)
    for i = 1:n
        if length(model[i].reactants) > 5 || length(model[i].products) > 5
            error("max: 5 reactants and 5 products ")
        end
    end

    temp_species = Array{Symbol,1}()
    temp_model = deepcopy(model)

    foreach(i -> append!(temp_species, append!(i.reactants, i.products)), temp_model)

    species = unique(temp_species)

    if :NULL == species[1]
        return species
    else
        species = [:NULL;species]
        species = unique(species)
    end
    return species
end

function getStoichio(model)
    n = length(model)
    stoichio = zeros(Int64, n, 11)
    for i = 1:n
        stoichio[i,1] = i
        if :NULL in model[i].reactants
            stoichio[i,2] = 0
            stoichio[i,3] = 0
            stoichio[i,4] = 0
            stoichio[i,5] = 0
            stoichio[i,6] = 0
        else
            if length(model[i].reactants) >= 1
                for j = 1:length(model[i].reactants)
                    if model[i].reactants[j] in model[i].products
                        ind = findfirst(x -> x == model[i].reactants[j],model[i].products)
                        if model[i].coeff_rea[j] > model[i].coeff_pro[ind]
                            stoichio[i,j+1] = - (model[i].coeff_rea[j] - model[i].coeff_pro[ind])
                        elseif model[i].coeff_rea[j] <= model[i].coeff_pro[ind]
                            stoichio[i,j+1] = 0
                        end
                    else
                        stoichio[i,j+1] = - model[i].coeff_rea[j]
                    end
                end
            end
        end
        if :NULL in model[i].products
            stoichio[i,7] = 0
            stoichio[i,8] = 0
            stoichio[i,9] = 0
            stoichio[i,10] = 0
            stoichio[i,11] = 0
        else
            if length(model[i].products) >= 1
                for j = 1:length(model[i].products)
                    if model[i].products[j] in model[i].reactants
                        ind = findfirst(x -> x == model[i].products[j],model[i].reactants)
                        if model[i].coeff_pro[j] > model[i].coeff_rea[ind]
                            stoichio[i,6+j] = model[i].coeff_pro[j] - model[i].coeff_rea[ind]
                        elseif model[i].coeff_pro[j] <= model[i].coeff_rea[ind]
                            stoichio[i,6+j] = 0
                        end
                    else
                        stoichio[i,6+j] =  model[i].coeff_pro[j]
                    end
                end
            end
        end
    end
    return stoichio
end

function getStoichio2(model)
    n = length(model)
    stoichio2 = zeros(Int64, n, 11)
    for i = 1:n
        stoichio2[i,1] = i
        if :NULL in model[i].reactants
            stoichio2[i,2] = 0
            stoichio2[i,3] = 0
            stoichio2[i,4] = 0
            stoichio2[i,5] = 0
            stoichio2[i,6] = 0
        elseif length(model[i].reactants) == 1
            stoichio2[i,2] = - model[i].coeff_rea[1]
        elseif length(model[i].reactants) == 2
            stoichio2[i,2] = - model[i].coeff_rea[1]
            stoichio2[i,3] = - model[i].coeff_rea[2]
        elseif length(model[i].reactants) == 3
            stoichio2[i,2] = - model[i].coeff_rea[1]
            stoichio2[i,3] = - model[i].coeff_rea[2]
            stoichio2[i,4] = - model[i].coeff_rea[3]
        elseif length(model[i].reactants) == 4
            stoichio2[i,2] = - model[i].coeff_rea[1]
            stoichio2[i,3] = - model[i].coeff_rea[2]
            stoichio2[i,4] = - model[i].coeff_rea[3]
            stoichio2[i,5] = - model[i].coeff_rea[4]
        elseif length(model[i].reactants) == 5
            stoichio2[i,2] = - model[i].coeff_rea[1]
            stoichio2[i,3] = - model[i].coeff_rea[2]
            stoichio2[i,4] = - model[i].coeff_rea[3]
            stoichio2[i,5] = - model[i].coeff_rea[4]
            stoichio2[i,6] = - model[i].coeff_rea[5]
        end

        if :NULL in model[i].products
            stoichio2[i,7] = 0
            stoichio2[i,8] = 0
            stoichio2[i,9] = 0
            stoichio2[i,10] = 0
            stoichio2[i,11] = 0
        elseif length(model[i].products) == 1
            stoichio2[i,7] = model[i].coeff_pro[1]
        elseif length(model[i].products) == 2
            stoichio2[i,7] = model[i].coeff_pro[1]
            stoichio2[i,8] = model[i].coeff_pro[2]
        elseif length(model[i].products) == 3
            stoichio2[i,7] = model[i].coeff_pro[1]
            stoichio2[i,8] = model[i].coeff_pro[2]
            stoichio2[i,9] = model[i].coeff_pro[3]
        elseif length(model[i].products) == 4
            stoichio2[i,7] = model[i].coeff_pro[1]
            stoichio2[i,8] = model[i].coeff_pro[2]
            stoichio2[i,9] = model[i].coeff_pro[3]
            stoichio2[i,10] = model[i].coeff_pro[4]
        elseif length(model[i].products) == 5
            stoichio2[i,7] = model[i].coeff_pro[1]
            stoichio2[i,8] = model[i].coeff_pro[2]
            stoichio2[i,9] = model[i].coeff_pro[3]
            stoichio2[i,10] = model[i].coeff_pro[4]
            stoichio2[i,11] = model[i].coeff_pro[5]
        end
    end
    return stoichio2
end


function getConnectivity(model,species)
    n = length(model)
    connectivity = zeros(Int64, n, 11)
    for i = 1:n
        connectivity[i,1] = i
        if length(model[i].reactants) == 1
            connectivity[i,2] = findfirst(x -> x == model[i].reactants[1], species)
        elseif length(model[i].reactants) == 2
            connectivity[i,2] = findfirst(x -> x == model[i].reactants[1], species)
            connectivity[i,3] = findfirst(x -> x == model[i].reactants[2], species)
        elseif length(model[i].reactants) == 3
            connectivity[i,2] = findfirst(x -> x == model[i].reactants[1], species)
            connectivity[i,3] = findfirst(x -> x == model[i].reactants[2], species)
            connectivity[i,4] = findfirst(x -> x == model[i].reactants[3], species)
        elseif length(model[i].reactants) == 4
            connectivity[i,2] = findfirst(x -> x == model[i].reactants[1], species)
            connectivity[i,3] = findfirst(x -> x == model[i].reactants[2], species)
            connectivity[i,4] = findfirst(x -> x == model[i].reactants[3], species)
            connectivity[i,5] = findfirst(x -> x == model[i].reactants[4], species)
        elseif length(model[i].reactants) == 5
            connectivity[i,2] = findfirst(x -> x == model[i].reactants[1], species)
            connectivity[i,3] = findfirst(x -> x == model[i].reactants[2], species)
            connectivity[i,4] = findfirst(x -> x == model[i].reactants[3], species)
            connectivity[i,5] = findfirst(x -> x == model[i].reactants[4], species)
            connectivity[i,6] = findfirst(x -> x == model[i].reactants[5], species)
        end

        if length(model[i].products) == 1
            connectivity[i,7] = findfirst(x -> x == model[i].products[1], species)
        elseif length(model[i].products) == 2
            connectivity[i,7] = findfirst(x -> x == model[i].products[1], species)
            connectivity[i,8] = findfirst(x -> x == model[i].products[2], species)
        elseif length(model[i].products) == 3
            connectivity[i,7] = findfirst(x -> x == model[i].products[1], species)
            connectivity[i,8] = findfirst(x -> x == model[i].products[2], species)
            connectivity[i,9] = findfirst(x -> x == model[i].products[3], species)
        elseif length(model[i].products) == 4
            connectivity[i,7] = findfirst(x -> x == model[i].products[1], species)
            connectivity[i,8] = findfirst(x -> x == model[i].products[2], species)
            connectivity[i,9] = findfirst(x -> x == model[i].products[3], species)
            connectivity[i,10] = findfirst(x -> x == model[i].products[4], species)
        elseif length(model[i].products) == 5
            connectivity[i,7] = findfirst(x -> x == model[i].products[1], species)
            connectivity[i,8] = findfirst(x -> x == model[i].products[2], species)
            connectivity[i,9] = findfirst(x -> x == model[i].products[3], species)
            connectivity[i,10] = findfirst(x -> x == model[i].products[4], species)
            connectivity[i,11] = findfirst(x -> x == model[i].products[5], species)
        end
    end
    return connectivity
end

function getInitiale(initiale,species)

    temp_initiale = deepcopy(initiale)
    ind = zeros(Int64,length(species),1)

    for i = 1:length(species)
        ind[i] = findfirst(x -> x == temp_initiale[i,1],species)
    end
    temp_initiale[:,1] = ind
    int_pop = temp_initiale
    return int_pop
end

function growth_estimate(growth_data)
    exp_growth(t, p) = p0 * exp.(p[1] * t)
    p0 = growth_data[1,5]
    fit = LsqFit.curve_fit(exp_growth, growth_data[:,1], growth_data[:,5], [0.0])
    param = fit.param[1]
    return param
end
