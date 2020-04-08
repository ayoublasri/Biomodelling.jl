function network_inference(genes::Array{String,1},path::String,sep::Char)

    sc_data = readdlm(path,sep)
    outt=bayNormJL.bayNorm(Data=sc_data[2:end,2:end],BETA_vec=nothing,Conditions = nothing,UMI_sffl = nothing,Prior_type = nothing,mode_version =false,mean_version=true,S = 20,FIX_MU = true,BB_SIZE_par = true, verbose = true)
    sc_data[2:end,2:end] = outt["Bay_out"]

    list_of_genes = convert(Array{String,1},sc_data[:,1])
    genes_index = zeros(Int,length(genes)+1)
    genes_index[1] = 1
    for i = 1:length(genes)
        genes_index[i+1] = findall(genes[i],list_of_genes, Levenshtein())[1]
    end

    if length(genes_index) != length(genes)+1
    return Error
    end

    sc_data = sc_data[genes_index,:]
    cd(pwd())
    writedlm("sc_data.csv",sc_data)

    number_of_genes = length(genes)
    dataset_name = "sc_data.csv"

    algorithm = PIDCNetworkInference()

    threshold = 0.15

    gene = get_nodes(dataset_name)

    network = InferredNetwork(algorithm, gene)

    adjacency_matrix1, labels_to_ids, ids_to_labels = get_adjacency_matrix(network, threshold)
    graph = LightGraphs.SimpleGraphs.SimpleGraph(adjacency_matrix1)

    number_of_nodes = size(adjacency_matrix1)[1]
    nodelabels = []
    for i in 1 : number_of_nodes
        push!(nodelabels, ids_to_labels[i])
    end
    return graph,nodelabels
end

# Example
# genes = ["MGMT";"FN1";"GLI1";"NUAK1";"TP53"]
# path = "-/Desktop/test_data.csv"
# sep = ','
# graphplot(graph,names=nodelabels,nodeshape=:circle)
