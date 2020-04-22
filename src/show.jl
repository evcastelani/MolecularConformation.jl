function Base.show(io::IO, c::ConformationOutput)
	print(io,"\n * Solver = $(c.solver)")
	print(io,"\n * Number of solutions = $(c.number)\n")
	print(io,  " * Solutions \n")	
	n = length(c.molecules[1].atoms)
	for k=1:c.number
		print(io,"   Molecule $(k) with LDE = $(c.molecules[k].lde) \n")
		for i=1:5
			print(io,"    ( $(c.molecules[k].atoms[i].x) , $(c.molecules[k].atoms[i].y) ,  $(c.molecules[k].atoms[i].z) ) \n")
		end
			print(io,"    ⋮  \n")
		for i=n-2:n
			print(io,"    ( $(c.molecules[k].atoms[i].x) , $(c.molecules[k].atoms[i].y) ,  $(c.molecules[k].atoms[i].z) ) \n")
		end
	
	end
	print(io," * Number of main operations [+-,*,/,√] = $(c.nop) \n")
	print(io," * Number of branchs  = $(c.nbranch) \n")
	print(io," * Number of pruning  = $(c.nprune) \n")
end
