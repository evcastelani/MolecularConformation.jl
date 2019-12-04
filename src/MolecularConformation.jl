module MolecularConformation
		
	export nmr,NMRinfo,NMRType,conformation,ConformationSetup,ConformationOutput,
			AtomType,MoleculeType,classical_bp, quaternion_bp,≈,
			generate_virtual_path
	
	
	# clean logging output
	using Logging
	conformation_logger = ConsoleLogger(stdout, Logging.Error)

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
	julia> options = ConformationSetup(0.001,5.5,classical_bp,true)

	julia> conformation(NMRdata,options)
	```
	as return a ConformationOutput type is provided.
	"""
	function conformation(NMRdata::NMRType,
			           cs::ConformationSetup)
		
#	print(" Cutting off distances greater than $(cs.cutoff) ... ")
#	dcutoff = findall(map(x->x>cs.cutoff,NMRdata.upperbound)) 
#	if isempty(dcutoff)
#		print("there are no distances greater than $(cs.cutoff) \n")
#	else
#		# Ainda precisamos discutir melhor a questão dos cortes
#		k = 0
#		for i in dcutoff
#			if NMRdata.label1 == "H" && NMRfile.label2 == "H"
#				NMRdata.lowerbound[i] = 0.0
#				NMRdata.upperbound[i] = 0.0
#				k += 1
#			end
#		end
#		print("$(k) distances were removed \n")
#	end
#	
#	print(" Checking if the file is a 3-click ... ")
#	n=last(NMRdata.vertex1)
#	virtual_path = generate_virtual_path(NMRdata)
#	if n == length(virtual_path)
#		print("the file is a 3-click \n")
#	else
#		print("using a virtual path for re-order \n")
#	end
#
#	
#
	print(" Solving the problem with $(cs.solver) ...")
#	
	solutions, t, bytes, gctime, memallocs = @timed cs.solver(NMRdata,cs.precision,cs.allsolutions)

	print(" Done! \n")
	end				

end
