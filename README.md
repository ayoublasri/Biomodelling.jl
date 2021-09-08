| **Build Status** | **Help** |
|:---:|:---:|
| [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] | [![][slack-img]][slack-url] |

# Biomodelling.jl

Authors:
- Ayoub Lasri (lasriay@gmail.com)
- Marc Sturrock

Framework for stochastic modelling in systems biology

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

Below outlines some typical way to define reactions

```julia 
    k1 = 1.0
    k2 = 0.2
    k3 = 10.0
    k4 = 0.1
    reaction1 = (name = "transcription", rate = k1, reactants = [:NULL], products =[:mRNA] , coeff_rea = [1] , coeff_pro = [1] )
    reaction2 = (name = "mRNA decay", rate = k2, reactants = [:mRNA], products =[:NULL], coeff_rea = [1], coeff_pro = [1])
    reaction3 = (name = "translation", rate = k3, reactants = [:mRNA], products =[:mRNA,:protein], coeff_rea = [1] , coeff_pro = [1,1] )
    reaction4 = (name = "protein decay", rate = k4, reactants = [:protein], products = [:NULL], coeff_rea = [1] , coeff_pro = [1] )
    model = (reaction1, reaction2, reaction3, reaction4)
```
Initial conditions and building the model

```julia 
    initiale_population = [:NULL 0;:A 10;:B 10]
    data = Biomodelling.Donne(model,initiale_population,1000,1.0, 1000, 0.38)
```
Model simulation: different algorithms were implemented, below an example using SSA.

```julia 
    T,X =Biomodelling.ssa(data)
```
