module MolecularConformation
		
	export conformation,ConformationSetup,ConformationOutput,
			AtomType,MoleculeType,classical_bp, quaternion_bp,≈
	
	
	# clean logging output
	using Logging
	conformation_logger = ConsoleLogger(stdout, Logging.Error)
	classical_bp_logger = ConsoleLogger(stdout, Logging.Error)
	quaternion_bp_logger = ConsoleLogger(stdout,Logging.Error)

	# loading basic packages
	using LinearAlgebra
	#using Quaternions
	import Base.show
	
	include("utils.jl")
	include("solvers.jl")
	include("show.jl")
	
	
	# main function to run the solver and related problem
	"""
	``` 
	conformation(D::Array{Float64,2},cs::ConformationSetup)
	
	```
	This is the main function in order to get the conformation of a molecule 
	from a distance matrix. To run this function, we need to setup a distance 
	array and then define some options using the ConformationSetup type. 

	## Example 
	
	```julia-repl
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
	julia> options = ConformationSetup(0.001,5.5,classical_bp,true)

	julia> conformation(D,options)
	```
	as return a ConformationOutput type is provided.
	"""
	function conformation(A::Array{Float64,2},cs::ConformationSetup)
		D=copy(A)
		print("\n Checking symmetry...")
		if D!=D'
			error("Distance matrix is not symmetric")
		else
			print(" Done!\n")
		end
		print(" Cut off distances greater than $(cs.cutoff)...")	
		(m,n) = size(D)
		nad = 0 #number of additional distances
		for i=1:n
			for j=i+4:n
				if D[i,j]>cs.cutoff
					D[i,j] = 0.0
					D[j,i] = 0.0
				else
					nad = nad + 1
				end
			end
		end
		print(" Done!\n")
		with_logger(conformation_logger) do
			@info "Distance Matrix" D
		end
		print(" Solving the problem with $(cs.solver) ...")
		
		solutions, t, bytes, gctime, memallocs = @timed cs.solver(n,D,nad,cs.precision,cs.allsolutions)
		print(" Done! \n")
		print(" Computing LDE for all solutions ...")
		worstlde = 0.0
		for i=1:solutions[1]
			
			solutions[2][i].lde = LDE(solutions[2][i],D,n,nad)
			if worstlde < solutions[2][i].lde
				worstlde = solutions[2][i].lde
			end
	#		println(storage_mol[i].atoms[n])
		end
		print("Done! \n")
		print(" Info: $(solutions[1]) solutions found with worst LDE given by $(worstlde) \n")
		return ConformationOutput(solutions[1],solutions[2],t,bytes,gctime)
		
	end

	function performance_conformation(A::Array{Float64,2},cs::ConformationSetup,ndiag::Int)
		D=copy(A)
		print("\n Checking symmetry...")
		if D!=D'
			error("Distance matrix is not symmetric")
		else
			print(" Done!\n")
		end
		print(" No cut off distances ...")	
		(m,n) = size(D)
		nad = 0 #number of additional distances
		for i=1:n
			for j=i+ndiag:n
					D[i,j] = 0.0
					D[j,i] = 0.0
		#			nad = nad + 1
			end
		end
		print(" Done!\n")
		with_logger(conformation_logger) do
			@info "Distance Matrix" D
		end
		print(" Solving the problem with $(cs.solver) ...")
		
		solutions, t, bytes, gctime, memallocs = @timed cs.solver(n,D,nad,cs.precision,cs.allsolutions)
		print(" Done! \n")
		print(" Computing LDE for all solutions ...")
		worstlde = 0.0
		for i=1:solutions[1]
			
			solutions[2][i].lde = LDE(solutions[2][i],D,n,nad)
			if worstlde < solutions[2][i].lde
				worstlde = solutions[2][i].lde
			end
	#		println(storage_mol[i].atoms[n])
		end
		print("Done! \n")
		print(" Info: $(solutions[1]) solutions found with worst LDE given by $(worstlde) \n")
		return ConformationOutput(solutions[1],solutions[2],t,bytes,gctime)
		
	end
end
