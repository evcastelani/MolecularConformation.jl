module MolecularConformation

export  preprocessing,NMRInfo,NMRType,conformation,ConformationSetup,ConformationOutput,
AtomType,MoleculeType,classicBP, classicBP_closure,
quaternionBP,symBP,≈,generate_virtual_path, bondangle,torsionmatrix,badtorsionangle,
pruningtest,LDE,build_distance_matrix,outputfilter,writefile,convert_to_dataframe


# loading basic packages
using LinearAlgebra,DelimitedFiles,SparseArrays,DataFrames,Dates
#using Quaternions
import Base.show

include("utils.jl")
include("solvers.jl")
include("show.jl")


# main function to run the solver and related problem
"""
``` 
conformation(NMRdata,cs::ConformationSetup)

```
This is the main function in order to get the conformation of a molecule. To run this function, we need to setup a NMR file and then define some options using the ConformationSetup type. 

## Example 

```julia-repl
julia> options = ConformationSetup(0.001,classicBP,true)

julia> conformation(NMRdata,options)
```
as return a ConformationOutput type is provided.

There are others parameters to setup, for example, need to complete.
"""
function conformation(NMRdata::NMRType,
		cs::ConformationSetup, 
		args...; 
		kargs...)


	s = cs.solver(NMRdata,cs.precision,cs.virtual_precision,cs.allsolutions, args...; kargs...)
	if cs.evalLDE == true
		map(i->MolecularConformation.LDE((typeof(s) == ConformationOutput ? s : s[1]).molecules[i],NMRdata),[1:1:(typeof(s) == ConformationOutput ? s : s[1]).number;])
	end
	return s
end				

end
