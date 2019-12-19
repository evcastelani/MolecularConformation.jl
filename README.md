# MolecularConformation.jl

![Build Status](https://travis-ci.com/evcastelani/MolecularConformation.jl.svg?branch=newinput)
![codecov.io](http://codecov.io/github/evcastelani/MolecularConformation.jl/coverage.svg?branch=newinput)

`MolecularConformation.jl` is a Julia Package which can be used to solve the conformation problem of some structure defined by a distance matrix of exact distance.

## Install 


This package was developed to Julia version >=1.0. Consequently, in order to install this package we need to type "]" in Julia REPL and after 

```julia
Pkg> add https://github.com/evcastelani/MolecularConformation.jl#nameofbranch
```

## Basic usage


This new version of Molecular Conformation works using a data list in `NMRType` format. In this sense, in order to run the `conformation function` we need to read a `.pdb` file. For example, let us consider the `pdb1a03.pdb`. In this case, we need convert the `.pdb` file to NMRType as follows 
Considering D a distance matrix given by

```julia
julia> data = nmr("pdb1a03.pdb") 
    
```
  
 Now,  we need to define the options of MolecularConformation.jl as 
  
```julia
julia> options = ConformationSetup(0.000001,classicBP,true)
```
where `0.000001` is the precision,  `classicBP` is the solver used to solve the problem and `true` value is used to compute all possible solutions.  In order to determine the positions  of atoms we to run the main function:
 
```julia
julia> conformation(data,options)
```
As return a ConformationOutput type is provided with all required information.

## TODO

1. [] Include signs vector;
1. [] Include in the data list the information about both atoms not just one 
1. [] Include in the data list the symmetry in the information too;
1. [] Discuss preprocessing;
1. [] Optimize redundant plan;
1. [] Implement chirality;
1. [] Optimize torsion angle in the repeating process. 



