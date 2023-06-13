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


## TO DO

1. [X] Discuss about the standard extension `.pdb` or `.mdjeep` or another; 
2. [X] Include signs vector (this idea doesn't work);
3. [X] Include in the data list the information about both atoms not just one 
4. [X] Include in the data list the symmetry in the information too;
5. [X] Discuss preprocessing;
6. [ ] Optimize redundant plan;
7. [ ] Implement chirality;
8. [X] Optimize torsion angle in the repeating process;
9.  [X] Modify the quaternion version to new input and improvements;
10. [X] Run all examples and create a table of tests;
11. [X] Include a decent documentation using `Documenter.jl`
12. [X] Including in preprocessament step the diagonal line (i,i) distance (to remove unnecessary if in main loop) 
13. [ ] Clear useless functions of COAP2 branch
14. [ ] Check and rewrite if necessary some docstrings
15. [ ] Merge COAP2 branch to master
16. [ ] Tag the last version to v1.0
17. [ ] Check the unitary tests
18. [ ] Check the tutorial related to last version
19. [ ] Update compiled results (pdf + images)
20. [ ] Write some readme files related to example folder and subfolders


