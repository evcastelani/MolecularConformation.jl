function Base.show(io::IO, c::ConformationOutput)
	print(io,"\n * Number of solutions = $(c.number)\n")
	print(io,  " * Solutions \n")	
	n = length(c.molecules[1].atoms)
	for k=1:c.number
		print(io,"   Molecule $(k) with LDE = $(c.molecules[k].lde) \n")
		for i=1:n
			print(io,"$(c.molecules[k].atoms[i])) \n")
			print(io,"    ( $(c.molecules[k].atoms[i].x) , $(c.molecules[k].atoms[i].y) ,  $(c.molecules[k].atoms[i].z) ) \n")
			print(io,"      ( $(c.molecules[k].atoms[i].x̂.x) , $(c.molecules[k].atoms[i].x̂.y) ,  $(c.molecules[k].atoms[i].x̂.z) ) \n")
			print(io,"      ( $(c.molecules[k].atoms[i].ŷ.x) , $(c.molecules[k].atoms[i].ŷ.y) ,  $(c.molecules[k].atoms[i].ŷ.z) ) \n")
			print(io,"      ( $(c.molecules[k].atoms[i].ẑ.x) , $(c.molecules[k].atoms[i].ẑ.y) ,  $(c.molecules[k].atoms[i].ẑ.z) ) \n")
		end
	end
	print(io," * Elapsed Time = $(c.elapsedtime)\n")
	print(io," * Allocated Memory = $(c.bytes)\n")
	print(io," * Garbage Collector Time = $(c.gctime)\n")
	print(io," * Number of main operations = $(c.nop \n)")
end
