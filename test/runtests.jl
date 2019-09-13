using Test, Statistics, ExcelReaders, Biomodelling
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

    initiale_population = [:NULL 0;:mRNA 5;:protein 100]
    data = Biomodelling.Donne(model,initiale_population,100,0.1, 100, 0.03,0.03)
    return data.M
end

@test test_1() == 4

function test_2()

    k1 = 4e-5
    k2 = 50
    k3 = 10
    k4 = 25

    Reaction1 = (name = "trans", rate = k1, reactants = [:A ; :B], products =[:A] , coeff_rea = [2;1] , coeff_pro = [3] )
    Reaction2 = (name = "birth", rate = k2, reactants = [:NULL], products =[:A] , coeff_rea = [1] , coeff_pro = [1] )
    Reaction3 = (name = "death", rate = k3, reactants = [:A], products =[:NULL] , coeff_rea = [1] , coeff_pro = [1] )
    Reaction4 = (name = "transf", rate = k4, reactants = [:NULL], products =[:B] , coeff_rea = [1] , coeff_pro = [1] )
    model = (Reaction1, Reaction2, Reaction3, Reaction4)

    initiale_population = [:NULL 0;:A 10;:B 10]

    data = Biomodelling.Donne(model,initiale_population,1000,1.0, 1000, 0.38)
    return data.stoichio
end

@test test_2() == [0 -1 1 0; 0 0 1 0; -1 0 0 0; 0 0 1 0]

function test_3(b,d)

    init = Int(round(b/d))
    Reaction1 = (name = "birth", rate = b, reactants = [:NULL], products =[:MDR1] , coeff_rea = [1] , coeff_pro = [1] )
    Reaction2 = (name = "death", rate = d, reactants = [:MDR1], products =[:NULL] , coeff_rea = [1] , coeff_pro = [1] )

    model = (Reaction1, Reaction2)

    initiale_population = [:NULL 0;:MDR1 init]
    maxtime = 100.0
    ts = 1.0
    data = Biomodelling.Donne(model,initiale_population,maxtime,ts,1000,0.38)
    sol_switch = zeros(1000)
    for i = 1:1000
        C,D = Biomodelling.tauleapswitch(data,100)
        sol_switch[i] = mean(D[:,2])
    end
    return sol_switch
end
@test mean(test_3(100.0,1.0)) - 100.0 < 0.1

function test_4()
    k1 = 4e-5
    k2 = 50
    k3 = 10
    k4 = 25

    Reaction1 = (name = "trans", rate = k1, reactants = [:A ; :B], products =[:A] , coeff_rea = [2;1] , coeff_pro = [3] )
    Reaction2 = (name = "birth", rate = k2, reactants = [:NULL], products =[:A] , coeff_rea = [1] , coeff_pro = [1] )
    Reaction3 = (name = "death", rate = k3, reactants = [:A], products =[:NULL] , coeff_rea = [1] , coeff_pro = [1] )
    Reaction4 = (name = "transf", rate = k4, reactants = [:NULL], products =[:B] , coeff_rea = [1] , coeff_pro = [1] )
    model = (Reaction1, Reaction2, Reaction3, Reaction4)

    initiale_population = [:NULL 0;:A 10;:B 10]

    growth = ExcelReaders.readxlsheet("test/Book1.xlsx","Sheet1",nrows=7,ncols=5)
    data = Biomodelling.Donne(model,initiale_population,1000,1.0, 1000, growth,0.03)
    return data.growth_rate
end

@test test_4() == 0.3970472060379792

function test_5()
    k1 = 4e-5
    k2 = 50
    k3 = 10
    k4 = 25

    Reaction1 = (name = "trans", rate = k1, reactants = [:A ; :B], products =[:A] , coeff_rea = [2;1] , coeff_pro = [3] )
    Reaction2 = (name = "birth", rate = k2, reactants = [:NULL], products =[:A] , coeff_rea = [1] , coeff_pro = [1] )
    Reaction3 = (name = "death", rate = k3, reactants = [:A], products =[:NULL] , coeff_rea = [1] , coeff_pro = [1] )
    Reaction4 = (name = "transf", rate = k4, reactants = [:NULL], products =[:B] , coeff_rea = [1] , coeff_pro = [1] )
    model = (Reaction1, Reaction2, Reaction3, Reaction4)

    initiale_population = [:NULL 0;:A 10;:B 10]
    growth = 0.11
    data = Biomodelling.Donne(model,initiale_population,1000,0.1, 10, growth,0.03)
    return data.growth_rate
end

@test test_5() == 0.11
