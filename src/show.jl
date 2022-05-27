function Base.show(io::IO,ct::Counter)
	print(io,"\n    ↳ In nodes           =  $(ct.node)")
	print(io,"\n    ↳ In virtual path    =  $(ct.virtual_path) ")
	print(io,"\n    ↳ In DDF             =  $(ct.ddf)   ")
	print(io,"\n    ↳ Number of branches =  $(ct.branch) " )
	print(io,"\n    ↳ Number of pruning  =  $(ct.prune) ")
end

function Base.show(io::IO, c::ConformationOutput)
	print(io,"\n ✔ Solver = $(c.solver)")
	print(io,"\n ✔ Number of solutions = $(c.number)\n")
	if c.number <1
		return
	end
	print(io,  " ✔ Solutions \n")	
	n = length(c.molecules[1].atoms)
	for k=1:c.number
		print(io,"   ↳ Molecule $(k) with LDE = $(c.molecules[k].lde) \n")
		for i=1:5
			print(io,"    ( $(c.molecules[k].atoms[i].x) , $(c.molecules[k].atoms[i].y) ,  $(c.molecules[k].atoms[i].z) ) \n")
		end
			print(io,"    ⋮  \n")
		for i=n-2:n
			print(io,"    ( $(c.molecules[k].atoms[i].x) , $(c.molecules[k].atoms[i].y) ,  $(c.molecules[k].atoms[i].z) ) \n")
		end
	
	end
	print(io," ✔ Number of main operations [+-,*,/,√] ")
	print(io,"   $(c.nop) ")
end
