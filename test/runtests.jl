using Test, Biomodelling
@test 1==1
function test_1()
    k1 = 1.0
    k2 = .2
    k3 = 10.0
    k4 = 0.1

    reaction1 = (name = "transcription", rate = k1, reactants = [:NULL], products =[:mRNA] , coeff_rea = [1] , coeff_pro = [1] )
    reaction2 = (name = "mRNA decay", rate = k2, reactants = [:mRNA], products =[:NULL], coeff_rea = [1], coeff_pro = [1])
    reaction3 = (name = "translation", rate = k3, reactants = [:mRNA], products =[:mRNA,:protein], coeff_rea = [1] , coeff_pro = [1,1] )
    reaction4 = (name = "protein decay", rate = k4, reactants = [:protein], products = [:NULL], coeff_rea = [1] , coeff_pro = [1] )
    model = (reaction1, reaction2, reaction3, reaction4)

    initiale_population = [:NULL 0;:mRNA 5;:protein 500]
    data = Biomodelling.Donne(model,initiale_population,1000,1.0,0.03)
    return data.M
end

@test test_1() == 4
