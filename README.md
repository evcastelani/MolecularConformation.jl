# MolecularConformation.jl

![Build Status](https://travis-ci.com/evcastelani/MolecularConformation.jl.svg?branch=master)
![codecov.io](http://codecov.io/github/evcastelani/MolecularConformation.jl/coverage.svg?branch=master)

`MolecularConformation.jl` is a Julia Package which can be used to solve the conformation problem of some structure defined by a distance matrix of exact distance.

## Install 

This package was developed to Julia version >=1.0. Consequently, in order to install this package we need to type "]" in Julia REPL and after 

```julia
Pkg> add https://github.com/evcastelani/MolecularConformation.jl
```

## Basic usage

Considering D a distance matrix given by

```julia
julia> D = [0.0      1.526    2.49239  2.89143  3.48974  4.58574  4.37632  3.94086  2.60817  3.06475;
				1.526    0.0      1.526    2.49239  2.93393  4.35209  4.66009  4.41794  3.03012  3.63256;
				2.49239  1.526    0.0      1.526    2.49239  3.83961  4.34457  3.88905  2.53092  2.97774;
				2.89143  2.49239  1.526    0.0      1.526    2.49239  2.9232   2.50687  1.49017  2.48433;
				3.48974  2.93393  2.49239  1.526    0.0      1.526    2.49239  2.93393  2.48279  3.84346;
				4.58574  4.35209  3.83961  2.49239  1.526    0.0      1.526    2.49239  2.90194  4.30309;
				4.37632  4.66009  4.34457  2.9232   2.49239  1.526    0.0      1.526    2.49239  3.84092;
				3.94086  4.41794  3.88905  2.50687  2.93393  2.49239  1.526    0.0      1.526    2.49239;
				2.60817  3.03012  2.53092  1.49017  2.48279  2.90194  2.49239  1.526    0.0      1.526;  
				3.06475  3.63256  2.97774  2.48433  3.84346  4.30309  3.84092  2.49239  1.526    0.0;]      
```
  
  we need to define the options of MolecularConformation.jl as 
  
```julia
julia> options = ConformationSetup(0.001,5.5,classical_bp,true)
```
where `0.001` is the precision, `5.5` is the cut off measure, `classical_bp` is the solver used to solve the problem and `true` value is used to compute all possible solutions.  In order to determine the positions  of atoms we to run the main function:
 
```julia
julia> conformation(D,options)
```
As return a ConformationOutput type is provided with all required information.
