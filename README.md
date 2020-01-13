# MolecularConformation.jl

![Build Status](https://travis-ci.com/evcastelani/MolecularConformation.jl.svg?branch=master)
![codecov.io](http://codecov.io/github/evcastelani/MolecularConformation.jl/coverage.svg?branch=master)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://evcastelani.github.io/MolecularConformation.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://evcastelani.github.io/MolecularConformation.jl/dev)

`MolecularConformation.jl` is a Julia Package which can be used to solve the conformation problem of some structure defined by a distance matrix of exact distance.

## Install 




This package was developed to Julia version >=1.0. Consequently, in order to install this package we need to type "]" in Julia REPL and after 

```julia
Pkg> add https://github.com/evcastelani/MolecularConformation.jl#nameofbranch
```

## Basic usage





This new version of Molecular Conformation works using a data list in `NMRType` format. In this sense, in order to run the `conformation function` we need to read a `.pdb` or `.mdjeep` file. For example, let us consider the `pdb1a03.mdjeep`. In this case, we need convert the `.mdjeep` file to `NMRType` as follows 

```julia
julia> data = preprocessing("pdb1a03.mdjeep") 
    
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

## TO DO


1. [ ] Discuss about the standard extension `.pdb` or `.mdjeep` or another; 
1. [X] Include signs vector (this idea doesn't work);
1. [X] Include in the data list the information about both atoms not just one 
1. [X] Include in the data list the symmetry in the information too;
1. [ ] Discuss preprocessing;
1. [ ] Optimize redundant plan;
1. [ ] Implement chirality;
1. [ ] Optimize torsion angle in the repeating process;
1. [X] Modify the quaternion version to new input and improvements;
1. [ ] Run all examples and create a table of tests;
1. [X] Include a decent documentation using `Documenter.jl`.

