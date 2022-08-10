| **Build Status** | **Help** |
|:---:|:---:|
| [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] | [![][slack-img]][slack-url] |

# Biomodelling.jl

Authors:
- Ayoub Lasri (lasriay@gmail.com)
- Marc Sturrock

Framework for stochastic modelling in systems biology as published in [BMC Bioinformatics](https://link.springer.com/article/10.1186/s12859-022-04778-9)

Usage questions can be posted in:
[Julia Community](https://julialang.org/community/)

[slack-img]: https://img.shields.io/badge/chat-on%20slack-yellow.svg
[slack-url]: https://julialang.slack.com

[travis-img]: https://travis-ci.org/ayoublasri/Biomodelling.jl.svg?branch=master
[travis-url]: https://travis-ci.org/ayoublasri/Biomodelling.jl

[codecov-img]: https://codecov.io/gh/ayoublasri/Biomodelling.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/ayoublasri/Biomodelling.jl

# Installation

```julia 
add https://github.com/ayoublasri/Biomodelling.jl.git 
```

# Simple Usage

We describe how to simulate simple reactions, in this case:

- reaction1 describes transcription with a rate k1.
- reaction2 describes mRNA decay with a rate k2.
- reaction3 describes translation with a rate k3.
- reaction4 describes protein decay with a rate k4.

Each reaction is a Named Tuple that has six entries:

```julia
NamedTuple{(:name, :rate, :reactants, :products, :coeff_rea, :coeff_pro), Tuple{String, Int64, Vector{Symbol}, Vector{Symbol}, Vector{Int64}, Vector{Int64}}}
```
- name: the name of the reaction (String).
- rate: rate of the reaction (Int64).
- reactants: the reaction reactant(s) given as a vector of symbols.
- products: the reaction product(s) given as a vector of symbols.
- coeff_rea and coeff_pro refer to the coefficient of the reactants and the products respectively and given as vector of Int64.

# Simple example

```julia 
    k1 = 1.0
    k2 = 0.2
    k3 = 10.0
    k4 = 0.1
    reaction1 = (name = "transcription", rate = k1, reactants = [:NULL], products =[:mRNA] , coeff_rea = [1] , coeff_pro = [1] )
    reaction2 = (name = "mRNA decay", rate = k2, reactants = [:mRNA], products =[:NULL], coeff_rea = [1], coeff_pro = [1])
    reaction3 = (name = "translation", rate = k3, reactants = [:mRNA], products =[:mRNA,:protein], coeff_rea = [1] , coeff_pro = [1,1] )
    reaction4 = (name = "protein decay", rate = k4, reactants = [:protein], products = [:NULL], coeff_rea = [1] , coeff_pro = [1] )
```

The model is then definied as the combiantion of the reactions defined above

```julia 
    model = (reaction1, reaction2, reaction3, reaction4)
```

# Initial conditions

Next step is to define the initial conditions, for this example we consider that at t=0, there are 5 mRNA molecules and 100 protein molecules. The initiale conditions should be given as Array{Any, 2} with the rows corresponding to model species including the NULL.

```julia 
    initiale_population = [:NULL 0;:mRNA 5;:protein 100]
```

# Building the model

The function Donne is then used to build a structure that contains information about the model and can be updated in time. The function Donne takes as inputs: the model, initial conditions, the simulation time (1000.0 in this example), time step (0.1 in this example), number of cells (100 in this example) and cells' growth rate (0.03 in this example).

```julia 
    data = Biomodelling.Donne(model,initiale_population,1000.0,0.1,100,0.03)
```
# Model simulation

Different algorithms were implemented, below an example using SSA. T is a vector of time steps. The species are sorted in the same order as the initial_population, for example to access the number of mRNA molecules for the entire simulation time: X[:,2].

```julia 
    T,X =Biomodelling.ssa(data)
```
