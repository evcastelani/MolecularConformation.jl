# MolecularConformation.jl

[![Build Status](https://travis-ci.com/evcastelani/MolecularConformation.jl.svg?branch=master)](https://travis-ci.com/evcastelani/MolecularConformation.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/rxequ5lhhisgw196/branch/master?svg=true)](https://ci.appveyor.com/project/evcastelani/molecularconformation-jl/branch/master)
[![codecov](https://codecov.io/gh/evcastelani/MolecularConformation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/evcastelani/MolecularConformation.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://evcastelani.github.io/MolecularConformation.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://evcastelani.github.io/MolecularConformation.jl/dev)

`MolecularConformation.jl` is a Julia Package which can be used to solve the conformation problem of some structure defined by a distance matrix of exact distance.

## Install 




This package was developed to Julia version >=1.0. Consequently, in order to install this package we need to type "]" in Julia REPL and after 

```julia
Pkg> add https://github.com/evcastelani/MolecularConformation.jl#nameofbranch
```

## Basic usage






This new version of Molecular Conformation works using a data list in `NMRType` format. In this sense, in order to run the `conformation function` we need to read a `.nmr` file. For example, let us consider the `pdb1a03.nmr` given in `examples` folder. In this case, we need convert the `.nmr` file to `NMRType` as follows 

```julia
julia> data = preprocessing("pdb1a03.nmr") 
    
```
  
 Now,  we need to define the options of MolecularConformation.jl as 
  
```julia
julia> options = ConformationSetup(0.000001,classicBP,true)
```
where `0.000001` is the precision,  `classicBP` is the solver used to solve the problem and `true` value is used to compute all possible solutions.  In order to determine the positions  of atoms we can run the main function:
 
```julia
julia> conformation(data,options)
```
As return a `ConformationOutput` type is provided with all required information.

## Note

Massive tests were carried out to prove the potential of the implementations. If you are interested in reproducing the tests, it is highly recommended that you use the `performance` branch. In addition, for pre-processing of instances, it was used [HCProtCLI](https://github.com/caomem/HCProtCLI).


## Benchmarks

To execute the benchmarks, you need to have the `MolecularConformation` package correctly installed. Follow these steps:

### Real instances

1. Navigate to the `real_instances` folder (type `cd examples/real_instances/`).
1. Run `julia`.
1. Load the necessary script using `include("perform2.jl")`.
1. To run all selected instances, use the command `perform("w", list_of_problems=Array{String,1}(), ε=1.0e-5, time_limit=Second(60), benchmarkSeconds=4500, benchmarkSamples=2, improv = (c,q) -> q/c)`.
1. Finally, execute `perform("w", list_of_problems=["pdb2k2f","pdb2kbm","pdb2j0z","pdb2adl"], ε=1.0e-5, time_limit=Second(60), benchmarkSeconds=4500, minSamples=20000, improv = (c,q) -> q/c)` to perform the specific `list_of_problems`.

If you want to run `BP-All`, simply execute `performRMSD("w", list_of_problems = ["pdb2k2f", "pdb2kbm"], ε=1.0e-6, time_limit=Second(120), benchmarkSeconds=20000, minSamples=10000, improv= (m,q) -> q/m)`.

### Artificial instances

1. Go to `examples/lavor_instances/`
1. Load `include("performance.jl")`
1. runperf(benchmarkSamples=100000,benchmarkSeconds=300,ε=maxintfloat(),virtual_ε=maxintfloat())

## TO DO

1. [X] Discuss about the standard extension `.pdb` or `.mdjeep` or another; 
1. [X] Include signs vector (this idea doesn't work);
1. [X] Include in the data list the information about both atoms not just one 
1. [X] Include in the data list the symmetry in the information too;
1. [X] Discuss preprocessing;
1. [ ] Optimize redundant plan;
1. [ ] Implement chirality;
1. [ ] Optimize torsion angle in the repeating process;
1. [X] Modify the quaternion version to new input and improvements;
1. [X] Run all examples and create a table of tests;
1. [X] Include a decent documentation using `Documenter.jl`.
1. [X] Merge branch master with COAP
1. [ ] Revise the Documentation (and the readme in benchmarks folders)
1. [ ] Add in README informations about the consolidate benchmark files and a brief description about the benchmark folders
1. [ ] Update the package to last LTS Julia version
