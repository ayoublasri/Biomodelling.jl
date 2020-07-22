function random_network(TS_num::Int64,gene_num::Int64,cell_num::Int64,activation_num::Int64,inhibition_num::Int64,birth_rate::Float64,decay_rate::Float64,activation_rate::Float64,inhibition_rate::Float64,div_noise::Float64)

        global k1 = birth_rate
        global k2 = decay_rate
        global k3 = activation_rate
        global k4 = inhibition_rate

        global ii = 0
        while ii < gene_num
                global ii = ii + 1
                @eval $(Symbol("reaction$ii")) = (name = "birth$ii", rate = k1, reactants = [:NULL], products =[Symbol(:G,ii)] , coeff_rea = [1] , coeff_pro = [1] )
                @eval $(Symbol("reaction$(ii+gene_num)")) = (name = "decay$ii", rate = k2, reactants = [Symbol(:G,ii)], products =[:NULL] , coeff_rea = [1] , coeff_pro = [1] )
        end

        gene_list  = Symbol.(:G,(1:gene_num))
        global kk = 0
        while kk < activation_num
                global kk = kk + 1
                global genes = sample(gene_list,2,replace=false)
                @eval $(Symbol("reaction$(kk+2*gene_num)")) = (name = "activation$kk", rate = k3, reactants = [genes[1]], products =[genes[2]], coeff_rea = [1] , coeff_pro = [1] )
        end

        global tt = 0
        while tt < inhibition_num
                global tt = tt + 1
                global genes = sample(gene_list,2,replace=false)
                @eval $(Symbol("reaction$(tt+2*gene_num+activation_num)")) = (name = "inhibition$tt", rate = k4, reactants = [genes[1],genes[2]], products =[genes[1]], coeff_rea = [1,1] , coeff_pro = [1] )
        end

        model = []
        for ll = 1:2*gene_num+activation_num+inhibition_num
                push!(model,@eval $(Symbol("reaction$ll")))
        end

        for i = 1:TS_num
                gene_list  = Symbol.(:G,(1:gene_num))
                initiale_population = [:NULL 0;gene_list rand(0:100,gene_num)]
                data = Biomodelling.Donne(model,initiale_population,1.0,0.2, cell_num, 0.03,0.03)
                temps, V, X = Biomodelling.exponential_growth(data,div_noise,ssa)

                table =  Array{Any,2}(undef,length(temps)+1,gene_num+1)
                table[1,1] = "time"
                table[2:end,1] = temps
                for l = 1:gene_num
                        table[1,l+1] = "G$l"
                        table[2,l+1] = data.X[l+1] ./ max(maximum(data.X),maximum(X))
                        table[3:end,l+1] = mean(X[:,:,l+1],dims=2)[2:end] ./ max(maximum(data.X),maximum(X))
                end
                writedlm("time_serie_$i.csv",table)
        end

        initiale_population = [:NULL 0;gene_list rand(0:100,gene_num)]
        data = Biomodelling.Donne(model,initiale_population,1000.0,0.1, cell_num, 0.03,0.03)
        temps, V, X = Biomodelling.exponential_growth(data,div_noise,ssa)

        table1 =  Array{Any,2}(undef,gene_num+1,cell_num+1)
        table1[1,1] = "time"
        table1[2:end,1] = gene_list
        for l = 1:cell_num
                table1[1,l+1] = "Cell$l"
        end
        for t = 1:gene_num
                table1[t+1,2:end] = X[end,:,t+1]
        end
        writedlm("SS_data.csv",table1)
        return model
end
