module MolecularConformation
		
	export  preprocessing,NMRInfo,NMRType,conformation,ConformationSetup,ConformationOutput,
			AtomType,MoleculeType,classicBP, classicBP_closure,
			quaternionBP,≈,generate_virtual_path, bondangle,torsionangle,torsionmatrix,
			pruningtest,LDE,build_distance_matrix,outputfilter
	
	
	# loading basic packages
	using LinearAlgebra,DelimitedFiles,SparseArrays
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
	"""
	function conformation(NMRdata::NMRType,
			           cs::ConformationSetup)
		

	print(" Solving the problem with $(cs.solver) ... ")
	solutions, t, bytes, gctime, memallocs = @timed cs.solver(NMRdata,cs.precision,cs.virtual_precision,cs.allsolutions)
	print(" Done! \n")
	print(" Computing the LDE for all solutions ... ")
	s = ConformationOutput(solutions[1],solutions[2],t,bytes,gctime)
	map(i->MolecularConformation.LDE(s.molecules[i],NMRdata),[1:1:s.number;])
	print(" Done! \n")
	return s
	end				

end
