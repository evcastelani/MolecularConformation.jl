function Base.show(io::IO,ct::Counter)
	print(io,"\n    $(ct.node) <- In nodes")
	print(io,"\n    $(ct.virtual_path) <- In virtual path")
	print(io,"\n    $(ct.ddf) <- In Direct Distance Feasibility")
	print(io,"\n    $(ct.branch) <- Number of branches" )
	print(io,"\n    $(ct.prune) <- Number of pruning")
end

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
	print(io," * Number of main operations [+-,*,/,√] \n")
	print(io, "  $(c.nop) \n")
end
