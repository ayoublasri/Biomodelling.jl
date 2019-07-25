mutable struct Donne

    M::Int64
    N::Int64
    T::Int64
    tau::Float64
    start::Float64
    switch_steps::Int64
    epsilon::Float64
    NoJ::Int64
    X::Matrix{Int64}
    stoichio::Matrix{Int64}
    stoichio2::Matrix{Int64}
    kr::Matrix{Float64}
    cm_rea::Matrix{Int64}
    cm_spe::Matrix{Int64}
    hor::Matrix{Int64}
    mrr::Matrix{Int64}

#   rname: file name that describes the model.

# ---- OUTPUT PARAMETERS ----
#   data.M: number of reactions, [1x1]
#   data.N: number of species, [1x1]
#   data.X: initial population of the system(), [Nx1]
#   data.stoichio: stoichiometric matrix, [Mx4]
#   data.kr: rate constants, [Mx1]
#   data.delays: reaction delays [type of delay and uniform dirstibution parameters],[Mx3]
#   data.cm_rea: connectivity matrix for the reactions, [Mx4]
#   data.cm_spe: connectivity matrix for the species, [MxK+1]

    function Donne(model, initiale, T, tau, epsilon=0.03)

        M = length(model)
#####################################
#### load stoichoimetric matrix #####
#####################################

        N = length(getSpecies(model))

        tmp_stoichio = getStoichio(model)
        tmp_stoichio2 = getStoichio2(model)

        stoichio = zeros(M, 4)
        stoichio2 = zeros(M,4)

        for m = 1:M
            stoichio[tmp_stoichio[m,1], :] = tmp_stoichio[m, 2:5]
        end
        for m = 1:M
            stoichio2[tmp_stoichio2[m,1], :] = tmp_stoichio2[m, 2:5]
        end
#######################################
#### load reaction rate constants #####
#######################################

        kr = zeros(M, 1)
        for i = 1:M
            kr[i,1] = model[i].rate
        end

############################################
#### load reaction connectivity matrix #####
############################################
        species = getSpecies(model)
        tmp_cm_rea = getConnectivity(model,species)

        cm_rea = zeros(Int64, M, 4)
        for m = 1:M
            cm_rea[tmp_cm_rea[m,1], :] = tmp_cm_rea[m, 2:5]
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
        K = 10

        cm_spe = zeros(Int64, N, K)
        no_rea = zeros(Int64, N, 1)
        for m = 1:M
            for k = 1:4
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

        for j = 1:M # run for all reactions

            for i = 1:2 # for the reactant species

                spe_i = cm_rea[j, i] # species index

                max_reac_order = maximum(abs.(stoichio[j , 1:2])) # maximum order of reactants of reaction j
                reac_order = abs(stoichio[j , i]) # order of species i in reaction j

                if spe_i !== 0 && spe_i !== 1 # if i-th species is not zero or NULL

                    if hor[spe_i] < max_reac_order # compare the required reactants
                        hor[spe_i] = max_reac_order
                    end

                    if mrr[spe_i] < reac_order # compare the required reactants
                        mrr[spe_i] = reac_order
                    end
                end
            end
        end
        new(M, N, T, tau, start, switch_steps, epsilon, NoJ, X, stoichio, stoichio2,kr, cm_rea, cm_spe, hor, mrr)
    end
end


function getSpecies(model)
    n = length(model)
    for i = 1:n
        if length(model[i].reactants) > 2 || length(model[i].products) > 2
            error("max: 2 reactants and 2 products ")
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
    stoichio = zeros(Int64, n, 5)
    for i = 1:n
        stoichio[i,1] = i
        if :NULL in model[i].reactants
            stoichio[i,2] = 0
            stoichio[i,3] = 0
        else
            if model[i].reactants[1] in model[i].products
                ind = findfirst(x -> x == model[i].reactants[1],model[i].products)
                if model[i].coeff_rea[1] > model[i].coeff_pro[ind]
                    stoichio[i,2] = - (model[i].coeff_rea[1] - model[i].coeff_pro[ind])
                elseif model[i].coeff_rea[1] <= model[i].coeff_pro[ind]
                    stoichio[i,2] = 0
                end
            else
                stoichio[i,2] = - model[i].coeff_rea[1]
            end
            if length(model[i].reactants) > 1
                if model[i].reactants[2] in model[i].products
                    ind = findfirst(x -> x == model[i].reactants[2],model[i].products)
                    if model[i].coeff_rea[2] > model[i].coeff_pro[ind]
                        stoichio[i,3] = - (model[i].coeff_rea[2] - model[i].coeff_pro[ind])
                    elseif model[i].coeff_rea[2] <= model[i].coeff_pro[ind]
                        stoichio[i,3] = 0
                    end
                else
                    stoichio[i,3] = - model[i].coeff_rea[2]
                end
            else
                stoichio[i,3] = 0
            end
        end
        if :NULL in model[i].products
            stoichio[i,4] = 0
            stoichio[i,5] = 0
        else
            if model[i].products[1] in model[i].reactants
                ind = findfirst(x -> x == model[i].products[1],model[i].reactants)
                if model[i].coeff_pro[1] > model[i].coeff_rea[ind]
                    stoichio[i,4] = model[i].coeff_pro[1] - model[i].coeff_rea[ind]
                elseif model[i].coeff_pro[1] <= model[i].coeff_rea[ind]
                    stoichio[i,4] = 0
                end
            else
                stoichio[i,4] = model[i].coeff_pro[1]
            end
            if length(model[i].products) > 1
                if model[i].products[2] in model[i].reactants
                    ind = findfirst(x -> x == model[i].products[2],model[i].reactants)
                    if model[i].coeff_pro[2] > model[i].coeff_rea[ind]
                        stoichio[i,5] = model[i].coeff_pro[2] - model[i].coeff_rea[ind]
                    elseif model[i].coeff_pro[2] <= model[i].coeff_rea[ind]
                        stoichio[i,5] = 0
                    end
                else
                    stoichio[i,5] =  model[i].coeff_pro[2]
                end
            else
                stoichio[i,5] = 0
            end
        end
    end
    return stoichio
end

function getStoichio2(model)
    n = length(model)
    stoichio2 = zeros(Int64, n, 5)
    for i = 1:n
        stoichio2[i,1] = i
        if :NULL in model[i].reactants
            stoichio2[i,2] = 0
            stoichio2[i,3] = 0
        elseif length(model[i].reactants) == 1
            stoichio2[i,2] = - model[i].coeff_rea[1]
            stoichio2[i,3] = 0
        else
            stoichio2[i,2] = - model[i].coeff_rea[1]
            stoichio2[i,3] = - model[i].coeff_rea[2]
        end
        if :NULL in model[i].products
            stoichio2[i,4] = 0
            stoichio2[i,5] = 0
        elseif length(model[i].products) == 1
            stoichio2[i,4] = model[i].coeff_pro[1]
            stoichio2[i,5] = 0
        else
            stoichio2[i,4] = model[i].coeff_pro[1]
            stoichio2[i,5] = model[i].coeff_pro[2]
        end
    end
    return stoichio2
end


function getConnectivity(model,species)
    n = length(model)
    connectivity = zeros(Int64, n, 5)
    for i = 1:n
        connectivity[i,1] = i
        if length(model[i].reactants) == 1
            connectivity[i,2] = findfirst(x -> x == model[i].reactants[1], species)
            connectivity[i,3] = 0
        else
            connectivity[i,2] = findfirst(x -> x == model[i].reactants[1], species)
            connectivity[i,3] = findfirst(x -> x == model[i].reactants[2], species)
        end
        if length(model[i].products) == 1
            connectivity[i,4] = findfirst(x -> x == model[i].products[1], species)
            connectivity[i,5] = 0
        else
            connectivity[i,4] = findfirst(x -> x == model[i].products[1], species)
            connectivity[i,5] = findfirst(x -> x == model[i].products[2], species)
        end
    end
    return connectivity
end

function getInitiale(initiale,species)

    temp_initiale = initiale
    ind = zeros(Int64,length(species),1)

    for i = 1:length(species)
        ind[i] = findfirst(x -> x == temp_initiale[i,1],species)
    end
    temp_initiale[:,1] = ind
    int_pop = temp_initiale
return int_pop
end
