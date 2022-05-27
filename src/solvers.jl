function classicBPseq(NMRdata::NMRType,
		ε :: Float64,
		virtual_ε :: Float64,
		allmol :: Bool,maxl=25,bug=true)
	if allmol == true
		error("This solver is not prepared to find all solutions yet")
	end
	n = NMRdata.dim
	mol = MoleculeType(Vector{AtomType}(undef,n),0.0)
	for i=1:n
		mol.atoms[i] = AtomType(0.0,0.0,0.0)
	end
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()
	n_prune = 0
	n_branch = 0
	#count_nop = [+-,*,/,√]
	nop_node = [0,0,0,0]
	nop_ddf = [0,0,0,0]
	nop_vpath = [0,0,0,0]
	# first atom
	mol.atoms[1].element = NMRdata.info[1,:].nzval[1].atom1
	mol.atoms[1].x = 0.0
	mol.atoms[1].y = 0.0
	mol.atoms[1].z = 0.0
	#second atom
	mol.atoms[2].element = NMRdata.info[2,:].nzval[1].atom1
	mol.atoms[2].x = -NMRdata.info[1,2].dist
	mol.atoms[2].y = 0.0
	mol.atoms[2].z = 0.0
	# tird atom
	D12 = NMRdata.info[1,2].dist
	D13 = NMRdata.info[1,3].dist
	D23 = NMRdata.info[2,3].dist
	D14 = 0.0
	D24 = 0.0
	D34 = 0.0
	cθ,sθ = bondangle(D12,D13,D23)
	nop_node += [3,6,1,1]
	cω,sω = (0.0,0.0)
	mol.atoms[3].element = NMRdata.info[3,:].nzval[1].atom1
	mol.atoms[3].x = -D12+D23*cθ
	mol.atoms[3].y = D23*sθ
	mol.atoms[3].z = 0.0
	nop_node += [1,2,0,0]
	C = zeros(4,4)
	C[1,4] = mol.atoms[3].x
	C[2,4] = mol.atoms[3].y 
	C[3,4] = mol.atoms[3].z
	C[1,1] = cθ
	C[1,2] = sθ
	C[2,1] = sθ
	C[2,2] = -cθ
	C[3,3] = -1.0
	C[4,4] = 1.0

	l = 4 # branching starts at atom 4
	pos = 4 # position in virtual path
	explore_right_side = zeros(Bool,n)
	C_list = Array{Array{Float64,2}}(undef,n) # to acess level l it is need to put l-3
	#	println(C)
	C_list[3] = copy(C)
	#	C_before = zeros(4,4)
	B = Array{Array{Float64,2}}(undef,n) # to acess level l it is need to put l-3
	B[4] = zeros(4,4)
	first_ocor = -1*ones(Int,n)
	virtual_count = 0
	while l<=n && l>3
		#if 4<l <= maxl
		#   println("B matrix in level $(l-1) ")
		#   display(B[l-1])
		#   println("C matrix")
		#   display(C_list[l-1])
		#end

		#display(l)

		if l==maxl+1
			error("bla")
		end
		#		if first_ocor[l-1] == -1 && first_ocor[NMRdata.virtual_path[l-1]] == -1 
		#			first_ocor[l-1] = NMRdata.virtual_path[l-1]
		#		end
		# TODO: otimizar!
		pos = findfirst(x->x==l-1,NMRdata.virtual_path) +1
		#C_before = zeros(4,4)
		C_before = copy(C_list[l-1])
		keep = true
		while keep
			try
				D14 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos]].dist
			catch
				error("$([NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos]])")
				D14 = sqrt((mol.atoms[NMRdata.virtual_path[pos-3]].x - mol.atoms[NMRdata.virtual_path[pos]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-3]].y - mol.atoms[NMRdata.virtual_path[pos]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-3]].z - mol.atoms[NMRdata.virtual_path[pos]].z)^2)  
			end			
			try 
				D24 = NMRdata.info[NMRdata.virtual_path[pos-2],NMRdata.virtual_path[pos]].dist
			catch		
				D24 = sqrt((mol.atoms[NMRdata.virtual_path[pos-2]].x - mol.atoms[NMRdata.virtual_path[pos]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-2]].y - mol.atoms[NMRdata.virtual_path[pos]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-2]].z - mol.atoms[NMRdata.virtual_path[pos]].z)^2)   
			end			
			try
				D34 = NMRdata.info[NMRdata.virtual_path[pos-1],NMRdata.virtual_path[pos]].dist
			catch
				D34 = sqrt((mol.atoms[NMRdata.virtual_path[pos-1]].x - mol.atoms[NMRdata.virtual_path[pos]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-1]].y - mol.atoms[NMRdata.virtual_path[pos]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-1]].z - mol.atoms[NMRdata.virtual_path[pos]].z)^2)   
			end	
			try 
				D12 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos-2]].dist
			catch
				D12 = sqrt((mol.atoms[NMRdata.virtual_path[pos-3]].x - mol.atoms[NMRdata.virtual_path[pos-2]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-3]].y - mol.atoms[NMRdata.virtual_path[pos-2]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-3]].z - mol.atoms[NMRdata.virtual_path[pos-2]].z)^2)   
			end	
			try
				D13 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos-1]].dist
			catch
				D13 = sqrt((mol.atoms[NMRdata.virtual_path[pos-3]].x - mol.atoms[NMRdata.virtual_path[pos-1]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-3]].y - mol.atoms[NMRdata.virtual_path[pos-1]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-3]].z - mol.atoms[NMRdata.virtual_path[pos-1]].z)^2)   
			end	
			try 
				D23 = NMRdata.info[NMRdata.virtual_path[pos-2],NMRdata.virtual_path[pos-1]].dist
			catch
				D23 = sqrt((mol.atoms[NMRdata.virtual_path[pos-2]].x - mol.atoms[NMRdata.virtual_path[pos-1]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-2]].y - mol.atoms[NMRdata.virtual_path[pos-1]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-2]].z - mol.atoms[NMRdata.virtual_path[pos-1]].z)^2)   
			end	
			cθ,sθ = bondangle(D23,D24,D34)
			cω,sω = badtorsionangle(D12,D13,D14,D23,D24,D34)
			#			println("l value = $(l) and NMRdatavalue = $(NMRdata.virtual_path[pos]) in position $(pos)")
			#	@show cθ, sθ,cω,sω, D34


			B[l] = torsionmatrix(cθ,sθ,cω,sω,D34)
			# if explore_right_side[l] == true 
			# 	B[l] = torsionmatrix(B[l])
			# end
			if l==NMRdata.virtual_path[pos]
				nop_node += [0,7,0,0] #torsion matrix
				nop_node += [3,6,1,1] # bond angle
				nop_node += [10,20,4,2] # bad torsion angle
				C_list[l] = prodmatrix(C_before,B[l]) 
				nop_node += [24,33,0,0] 
				keep = false

			else
				nop_node += [0,7,0,0] #torsion matrix
				nop_node += [3,6,1,1] # bond angle
				nop_node += [10,20,4,2] # bad torsion angle

				nop_vpath += [0,7,0,0] # torsion matrix
				nop_vpath += [3,6,1,1] # bond angle
				nop_vpath += [10,20,4,2] # bad torsion angle

				cpx = mol.atoms[NMRdata.virtual_path[pos]].x
				cpy = mol.atoms[NMRdata.virtual_path[pos]].y
				cpz = mol.atoms[NMRdata.virtual_path[pos]].z
				if bug
					println("($cpx , $cpy ,$cpz )")
				end
				Virtual_Torsion = prodmatrix(C_before,B[l])
				# println("Virtual Torsion = $(C_before) * $(B[l])$(Virtual_Torsion)")
				nop_vpath += [24,33,0,0]
				nop_node += [23,33,0,0]

				if sqrt((Virtual_Torsion[1,4]- cpx)^2+(Virtual_Torsion[2,4]- cpy)^2+(Virtual_Torsion[3,4]- cpz)^2)> virtual_ε
					B[l] = torsionmatrix(B[l])
					C_before = prodmatrix(C_before,B[l])
					nop_vpath += [24,33,0,0]
					nop_node += [24,33,0,0]
					#	println("passou 1")
				else
					#println("passou 2")
					C_before = copy(Virtual_Torsion) 
					#C_before = Virtual_Torsion
				end
				if bug
					println("Virtual Torsion")
					display(Virtual_Torsion)
				end
				#println("Torsion matrix $(B[l])")
				#println("C_before matrix $(C_before)")
				#@debug "virtual atom position  " C_before[1,4],C_before[2,4],C_before[3,4]
				pos = pos+1		
			end
		end


		if explore_right_side[l] == false
			if bug && (3<l<= maxl) 
				println("B matrix in level $(l) in left side")
				display(B[l])
				println("C matrix ")
				display(C_list[l])
			end
			mol.atoms[l].element = NMRdata.info[l,:].nzval[1].atom1		
			mol.atoms[l].x = C_list[l][1,4]
			mol.atoms[l].y = C_list[l][2,4]
			mol.atoms[l].z = C_list[l][3,4]
			count = [0,0,0,0]
			λ , count  = pruningtest(mol,l,NMRdata,ε,count) 
			nop_ddf += count
			#println("C = C_before*B at level $(l) left side  $(C_list[l])) = $(C_list[l-1]) * $(B[l])")
			if λ == 1 
				if l<n
					#		println("Partial solution by left side at level $(l) " , mol)
					n_branch +=1
					#classicBP_closure(l+1,pos+1,mol,C)
					#explore_right_side[l+1]=false
				else
					nsol=nsol+1
					storage_mol[nsol] = copy(mol)
					@debug "Rank n was reached, a solution was found " 
				end
			else
				n_prune += 1
				explore_right_side[l] = true
			end
		end
		if explore_right_side[l] == true 
			#display(l)
			B[l] = torsionmatrix(B[l])
			#nop_node += [0,0,0,0]
			C_list[l] = prodmatrix(C_before,B[l])# tenho que otimizar este calculo
			if bug && (3<l<=maxl)
				println("B matrix in level $(l) in right side")
				display(B[l])
				println("C matrix")
				display(C_list[l])
			end
			nop_node += [24,33,0,0]
			mol.atoms[l].x = C_list[l][1,4]
			mol.atoms[l].y = C_list[l][2,4]
			mol.atoms[l].z = C_list[l][3,4]
			count = [0,0,0,0]
			ρ ,count = pruningtest(mol,l,NMRdata,ε,count) #preciso modificar
			nop_ddf += count 
			#println("C = C_before*B at level $(l) right side  $(C_list[l])) = $(C_list[l-1]) * $(B[l])")
			if ρ == 1 
				if l<n
					#		println("Partial solution by right side at level $(l) ", mol)
					n_branch += 1
					#explore_right_side[l+1]=false
					#classicBP_closure(l+1,pos+1,mol,C)
				else
					nsol = nsol+1
					storage_mol[nsol] = copy(mol)			
					@debug "Rank n was reached, a solution was found " 
				end
			else
				explore_right_side[l] = false
				k = l-1
				while explore_right_side[k] == true
					explore_right_side[k] = false
					k -= 1
				end
				explore_right_side[k] = true
				l = k-1
				n_prune += 1
			end


		end
		l += 1
	end
	if l == 3 
		error("Solution not found, problem possible infeasible")
	end
	display(first_ocor)
	return nsol, storage_mol,Counter(nop_node,nop_vpath,nop_ddf,n_branch,n_prune)
end

################
function classicBPseq2(NMRdata::NMRType,
		ε :: Float64,
		virtual_ε :: Float64,
		allmol :: Bool)
	if allmol == true
		error("This solver is not prepared to find all solutions yet")
	end

	n = NMRdata.dim
	mol = MoleculeType(Vector{AtomType}(undef,n),0.0)
	for i=1:n
		mol.atoms[i] = AtomType(0.0,0.0,0.0)
	end
	C = zeros(4,4)
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()
	n_prune = 0
	n_branch = 0
	#count_nop = [+-,*,/,√]
	nop_node = [0,0,0,0]
	nop_ddf = [0,0,0,0]
	nop_vpath = [0,0,0,0]
	#l = 1 first atom	# first atom
	mol.atoms[1].element = NMRdata.info[1,:].nzval[1].atom1
	mol.atoms[1].x = 0.0
	mol.atoms[1].y = 0.0
	mol.atoms[1].z = 0.0
	#second atom
	mol.atoms[2].element = NMRdata.info[2,:].nzval[1].atom1
	mol.atoms[2].x = -NMRdata.info[1,2].dist
	mol.atoms[2].y = 0.0
	mol.atoms[2].z = 0.0
	# tird atom
	D12 = NMRdata.info[1,2].dist
	D13 = NMRdata.info[1,3].dist
	D23 = NMRdata.info[2,3].dist
	D14 = 0.0
	D24 = 0.0
	D34 = 0.0
	cθ,sθ = bondangle(D12,D13,D23)
	nop_node += [3,6,1,1]
	cω,sω = (0.0,0.0)
	mol.atoms[3].element = NMRdata.info[3,:].nzval[1].atom1
	mol.atoms[3].x = -D12+D23*cθ
	mol.atoms[3].y = D23*sθ
	mol.atoms[3].z = 0.0
	nop_node += [1,2,0,0]
	C = zeros(4,4)
	C[1,4] = mol.atoms[3].x
	C[2,4] = mol.atoms[3].y 
	C[3,4] = mol.atoms[3].z
	C[1,1] = cθ
	C[1,2] = sθ
	C[2,1] = sθ
	C[2,2] = -cθ
	C[3,3] = -1.0
	C[4,4] = 1.0
	l = 4 # branching starts at atom 4
	pos = 4 # position in virtual path
	explore_right_side = zeros(Bool,n)
	C_list = Array{Array{Float64,2}}(undef,n) # to acess level l it is need to put l-3
	println(C)
	C_list[3] = C
	C_before = zeros(4,4)
	B = Array{Array{Float64,2}}(undef,n) # to acess level l it is need to put l-3
	B[4] = zeros(4,4)
	while l<=n
		pos = findall(NMRdata.virtual_path.==l-1)[1] +1
		display(l)
		copyto!(C_before,C_list[l-1])
		keep = true
		while keep
			try
				D14 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos]].dist
			catch
				D14 = sqrt((mol.atoms[NMRdata.virtual_path[pos-3]].x - mol.atoms[NMRdata.virtual_path[pos]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-3]].y - mol.atoms[NMRdata.virtual_path[pos]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-3]].z - mol.atoms[NMRdata.virtual_path[pos]].z)^2)  
			end			
			try 
				D24 = NMRdata.info[NMRdata.virtual_path[pos-2],NMRdata.virtual_path[pos]].dist
			catch		
				D24 = sqrt((mol.atoms[NMRdata.virtual_path[pos-2]].x - mol.atoms[NMRdata.virtual_path[pos]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-2]].y - mol.atoms[NMRdata.virtual_path[pos]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-2]].z - mol.atoms[NMRdata.virtual_path[pos]].z)^2)   
			end			
			try
				D34 = NMRdata.info[NMRdata.virtual_path[pos-1],NMRdata.virtual_path[pos]].dist
			catch
				D34 = sqrt((mol.atoms[NMRdata.virtual_path[pos-1]].x - mol.atoms[NMRdata.virtual_path[pos]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-1]].y - mol.atoms[NMRdata.virtual_path[pos]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-1]].z - mol.atoms[NMRdata.virtual_path[pos]].z)^2)   
			end	
			try 
				D12 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos-2]].dist
			catch
				D12 = sqrt((mol.atoms[NMRdata.virtual_path[pos-3]].x - mol.atoms[NMRdata.virtual_path[pos-2]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-3]].y - mol.atoms[NMRdata.virtual_path[pos-2]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-3]].z - mol.atoms[NMRdata.virtual_path[pos-2]].z)^2)   
			end	
			try
				D13 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos-1]].dist
			catch
				D13 = sqrt((mol.atoms[NMRdata.virtual_path[pos-3]].x - mol.atoms[NMRdata.virtual_path[pos-1]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-3]].y - mol.atoms[NMRdata.virtual_path[pos-1]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-3]].z - mol.atoms[NMRdata.virtual_path[pos-1]].z)^2)   
			end	
			try 
				D23 = NMRdata.info[NMRdata.virtual_path[pos-2],NMRdata.virtual_path[pos-1]].dist
			catch
				D23 = sqrt((mol.atoms[NMRdata.virtual_path[pos-2]].x - mol.atoms[NMRdata.virtual_path[pos-1]].x)^2 + (mol.atoms[NMRdata.virtual_path[pos-2]].y - mol.atoms[NMRdata.virtual_path[pos-1]].y)^2+ (mol.atoms[NMRdata.virtual_path[pos-2]].z - mol.atoms[NMRdata.virtual_path[pos-1]].z)^2)   
			end	
			cθ,sθ = bondangle(D23,D24,D34)
			cω,sω = badtorsionangle(D12,D13,D14,D23,D24,D34)
			println("l value = $(l) and NMRdatavalue = $(NMRdata.virtual_path[pos]) in position $(pos)")
			B[l] = torsionmatrix(cθ,sθ,cω,sω,D34)
			# if explore_right_side[l] == true 
			# 	B[l] = torsionmatrix(B[l])
			# end
			if l==NMRdata.virtual_path[pos]
				nop_node += [0,7,0,0] #torsion matrix
				nop_node += [3,6,1,1] # bond angle
				nop_node += [10,20,4,2] # bad torsion angle
				C_list[l] = prodmatrix(C_list[l-1],B[l]) 
				nop_node += [24,33,0,0] 
				keep = false
			else
				nop_node += [0,7,0,0] #torsion matrix
				nop_node += [3,6,1,1] # bond angle
				nop_node += [10,20,4,2] # bad torsion angle

				nop_vpath += [0,7,0,0] # torsion matrix
				nop_vpath += [3,6,1,1] # bond angle
				nop_vpath += [10,20,4,2] # bad torsion angle

				cpx = mol.atoms[NMRdata.virtual_path[pos]].x
				cpy = mol.atoms[NMRdata.virtual_path[pos]].y
				cpz = mol.atoms[NMRdata.virtual_path[pos]].z
				Virtual_Torsion = prodmatrix(C_list[l-1],B[l])
				println("Matriz Virtual Torsion")
				display(Virtual_Torsion)
				# println("Virtual Torsion = $(C_before) * $(B[l])$(Virtual_Torsion)")
				nop_vpath += [24,33,0,0]
				nop_node += [23,33,0,0]

				if sqrt((Virtual_Torsion[1,4]- cpx)^2+(Virtual_Torsion[2,4]- cpy)^2+(Virtual_Torsion[3,4]- cpz)^2)> virtual_ε
					B[l] = torsionmatrix(B[l])
					C_list[l-1] = prodmatrix(C_list[l-1],B[l])
					nop_vpath += [24,33,0,0]
					nop_node += [24,33,0,0]
					#println("passou 1")
				else
					#println("passou 2")
					copyto!(C_list[l-1],Virtual_Torsion) 
					#C_before = Virtual_Torsion
				end
				#println("Torsion matrix $(B[l])")
				#println("C_before matrix $(C_before)")
				#@debug "virtual atom position  " C_before[1,4],C_before[2,4],C_before[3,4]
				pos = pos+1		
			end
		end


		if explore_right_side[l] == true 
			B[l] = torsionmatrix(B[l])
			#nop_node += [0,0,0,0]
			C_list[l] = prodmatrix(C_list[l-1],B[l])# tenho que otimizar este calculo
			nop_node += [24,33,0,0]
		end
		mol.atoms[l].x = C_list[l][1,4]
		mol.atoms[l].y = C_list[l][2,4]
		mol.atoms[l].z = C_list[l][3,4]
		count = [0,0,0,0]
		ρ ,count = pruningtest(mol,l,NMRdata,ε,count) #preciso modificar
		nop_ddf += count 
		#println("C = C_before*B at level $(l) right side  $(C_list[l])) = $(C_list[l-1]) * $(B[l])")
		if ρ == 1 
			if l<n
				#		println("Partial solution by right side at level $(l) ", mol)
				n_branch += 1
				l += 1
			else
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)			
				@debug "Rank n was reached, a solution was found " 
			end
		else
			if explore_right_side[l] == true
				explore_right_side[l] = false
				k = l-1
				while explore_right_side[k] == true
					explore_right_side[k] = false
					k -= 1
				end
				explore_right_side[k] = true
				l = k
				n_prune += 1
			else
				explore_right_side[l] = true
			end
		end
	end
	return nsol, storage_mol,Counter(nop_node,nop_vpath,nop_ddf,n_branch,n_prune)
end


##############

"""
```
classicBP :: Function
```
This function is an important solver of this package. This solver implements the algorithm given in 

Liberti, L., Lavor, C., & Maculan, N. (2008). A branch‐and‐prune algorithm for the molecular distance geometry problem. International Transactions in Operational Research, 15(1), 1-17.

The main difference of our implementation is that the input data (`NMRType`) can store informations of preprocessing function and as consequence we can optimize the search tree of this algorithm. 
"""
function classicBP(NMRdata :: NMRType,
		ε :: Float64,
		virtual_ε :: Float64,
		allmol :: Bool, time_limit)
	
	start =  Dates.now()
	time_elapsed = Second(0.0)
	
	n = NMRdata.dim
	if n < 3
		ArgumentError("Invalid dimension of NMRdata")
	end

	virtual_ε² = virtual_ε*virtual_ε 
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()

	function initialization()
		mol = MoleculeType(Vector{AtomType}(undef,n),0.0)
		for i=1:n
			mol.atoms[i] = AtomType(0.0,0.0,0.0)
		end
		C = zeros(4,4)
			
		# first atom
		mol.atoms[1].element = NMRdata.info[1,:].nzval[1].atom1
		mol.atoms[1].x = 0.0
		mol.atoms[1].y = 0.0
		mol.atoms[1].z = 0.0
		#second atom
		mol.atoms[2].element = NMRdata.info[2,:].nzval[1].atom1
		mol.atoms[2].x = -NMRdata.info[1,2].dist
		mol.atoms[2].y = 0.0
		mol.atoms[2].z = 0.0
		# tird atom
		D12 = NMRdata.info[1,2].dist
		D13 = NMRdata.info[1,3].dist
		D23 = NMRdata.info[2,3].dist
		cθ,sθ = bondangle(D12,D13,D23)
		mol.atoms[3].element = NMRdata.info[3,:].nzval[1].atom1
		mol.atoms[3].x = -D12+D23*cθ
		mol.atoms[3].y = D23*sθ
		mol.atoms[3].z = 0.0
		C = zeros(4,4)
		C[1,4] = mol.atoms[3].x
		C[2,4] = mol.atoms[3].y 
		C[3,4] = mol.atoms[3].z
		C[1,1] = cθ
		C[1,2] = sθ
		C[2,1] = sθ
		C[2,2] = -cθ
		C[3,3] = -1.0
		C[4,4] = 1.0

		return 4,4,mol,C
	end
 
	# defining closure
	function classicBP_closure(l :: Int64,
		pos::Int64,
		mol :: MoleculeType,
		C :: Array{Float64,2})
		
		time_elapsed = Dates.now()-start
		if time_elapsed>time_limit && l<n
			error("Time limit reached without found a solution!")
		end

		C_before = copy(C)
		while true
			virtualPos = NMRdata.virtual_path[pos]
			virtualLastPos = NMRdata.virtual_path[pos-1]
			
			if NMRdata.virtual_path[pos-3] == virtualPos
				D14 = 0.0
			else
				D14 = NMRdata.info[NMRdata.virtual_path[pos-3],virtualPos].dist
			end
			if NMRdata.virtual_path[pos-2] == virtualPos
				D24 = 0.0
			else
				D24 = NMRdata.info[NMRdata.virtual_path[pos-2],virtualPos].dist
			end

			if virtualLastPos == virtualPos
				D34 = 0.0
			else
				D34 = NMRdata.info[virtualLastPos,virtualPos].dist
			end

			if NMRdata.virtual_path[pos-3] == NMRdata.virtual_path[pos-2]
				D12 = 0.0
			else
				D12 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos-2]].dist
			end
			if NMRdata.virtual_path[pos-3] == virtualLastPos
				D13 = 0.0
			else
				D13 = NMRdata.info[NMRdata.virtual_path[pos-3],virtualLastPos].dist
			end
			if NMRdata.virtual_path[pos-2] == virtualLastPos
				D23 = 0.0
			else
				D23 = NMRdata.info[NMRdata.virtual_path[pos-2],virtualLastPos].dist
			end

			cθ,sθ = bondangle(D23,D24,D34)
			cω,sω = badtorsionangle(D12,D13,D14,D23,D24,D34)
			B = torsionmatrix(cθ,sθ,cω,sω,D34)
			if l==virtualPos
				C = prodmatrix(C_before,B)
				break
			else
				cpx = mol.atoms[virtualPos].x
				cpy = mol.atoms[virtualPos].y
				cpz = mol.atoms[virtualPos].z

				Virtual_Torsion = prodmatrix(C_before,B)

				if (Virtual_Torsion[1,4]- cpx)^2+(Virtual_Torsion[2,4]- cpy)^2+(Virtual_Torsion[3,4]- cpz)^2> virtual_ε²
					B = torsionmatrix(B)
					C_before = prodmatrix(C_before,B)
				else
					C_before = copy(Virtual_Torsion) 
				end
				pos = pos+1		
			end
		end

		mol.atoms[l].element = NMRdata.info[l,:].nzval[1].atom1		
		mol.atoms[l].x = C[1,4]
		mol.atoms[l].y = C[2,4]
		mol.atoms[l].z = C[3,4]

		if pruningtest(mol,l,NMRdata,ε) 
			if l<n
				classicBP_closure(l+1,pos+1,mol,C)
			else
				nsol=nsol+1
				storage_mol[nsol] = copy(mol)
				return
			end
		end
		if !allmol && nsol>0
			return
		end

		B = torsionmatrix(B)
		C = prodmatrix(C_before,B)# tenho que otimizar este calculo
		mol.atoms[l].x = C[1,4]
		mol.atoms[l].y = C[2,4]
		mol.atoms[l].z = C[3,4]
		
		if pruningtest(mol,l,NMRdata,ε) #preciso modificar 
			if l<n
				classicBP_closure(l+1,pos+1,mol,C)
			else
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				return
			end
		end
	end # closure

	__l, __pos, __mol, __C = initialization()
	classicBP_closure(__l, __pos, __mol, __C)
	return nsol, storage_mol

end #solver classicBP

""" 
```
quaternionBP :: Function
```
This function defines a new solver. The implementation follows ideas describing in:

Fidalgo, F. Using Quaternion Geometric Algebra for efficient rotations in Branch and Prune Algorithtm to solve the Discretizable Molecular Distance Geometry Problem. In: Proceedings of AGACSE 2018, Campinas-SP, Brazil.
"""
function quaternionBP(NMRdata :: NMRType,
		ε :: Float64,
		virtual_ε :: Float64,
		allmol :: Bool, time_limit)
	
	start =  Dates.now()
	time_elapsed = Second(0.0)

	n = NMRdata.dim
	if n < 3
		ArgumentError("Invalid dimension of NMRdata")
	end

	virtual_ε² = virtual_ε*virtual_ε 
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()
	
	function initialization()
		# TODO: (Emerson) na criação desse vetor você não pode já estabelecer um valor default para atoms?
		mol = MoleculeType(Vector{AtomType}(undef,n),0.0)
		for i=1:n
			mol.atoms[i] = AtomType(0.0,0.0,0.0)
		end
		Q = Quaternion(0.0,0.0,0.0,0.0)

		# first atom
		mol.atoms[1].element = NMRdata.info[1,:].nzval[1].atom1
		mol.atoms[1].x = 0.0
		mol.atoms[1].y = 0.0
		mol.atoms[1].z = 0.0
		#second atom
		mol.atoms[2].element = NMRdata.info[2,:].nzval[1].atom1
		mol.atoms[2].x = -NMRdata.info[1,2].dist
		mol.atoms[2].y = 0.0
		mol.atoms[2].z = 0.0
		# tird atom
		D12 = NMRdata.info[1,2].dist
		D13 = NMRdata.info[1,3].dist
		D23 = NMRdata.info[2,3].dist
		cθ,sθ = qbondangle(D12,D13,D23)
		Q = Quaternion(0.0,-cθ,-sθ,0.0)
		d = 2.0*D23
		qmol = Quaternion(0.0,d*(cθ*cθ-0.5),d*(cθ*sθ),0.0)
		mol.atoms[3].element = NMRdata.info[3,:].nzval[1].atom1
		mol.atoms[3].x = qmol.v1 + mol.atoms[2].x
		mol.atoms[3].y = qmol.v2 + mol.atoms[2].y
		mol.atoms[3].z = qmol.v3 + mol.atoms[2].z

		return 4,4,mol,Q
	end
	
	# defining closure
	function quaternionBP_closure(l :: Int64,
									pos::Int64,
									mol :: MoleculeType,
									Q :: Quaternion)

		time_elapsed = Dates.now()-start
		if time_elapsed>time_limit && l<n
			error("Time limit reached without found a solution!")
		end
	
		lastpos = 1
		D34 = 0.0
		virtualLastPos = 0.0
		a = 0.0
		b = 0.0
		c = 0.0
		d = 0.0 
		Q_before = copy(Q)

		while true
			virtualPos = NMRdata.virtual_path[pos]
			virtualLastPos = NMRdata.virtual_path[pos-1]

			if NMRdata.virtual_path[pos-3] == virtualPos
				D14 = 0.0
			else
				D14 = NMRdata.info[NMRdata.virtual_path[pos-3],virtualPos].dist
			end
			if NMRdata.virtual_path[pos-2] == virtualPos
				D24 = 0.0
			else
				D24 = NMRdata.info[NMRdata.virtual_path[pos-2],virtualPos].dist
			end

			if virtualLastPos == virtualPos
				D34 = 0.0
			else
				D34 = NMRdata.info[virtualLastPos,virtualPos].dist
			end

			if NMRdata.virtual_path[pos-3] == NMRdata.virtual_path[pos-2]
				D12 = 0.0
			else
				D12 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos-2]].dist
			end
			if NMRdata.virtual_path[pos-3] == virtualLastPos
				D13 = 0.0
			else
				D13 = NMRdata.info[NMRdata.virtual_path[pos-3],virtualLastPos].dist
			end
			if NMRdata.virtual_path[pos-2] == virtualLastPos
				D23 = 0.0
			else
				D23 = NMRdata.info[NMRdata.virtual_path[pos-2],virtualLastPos].dist
			end
			
			cθ,sθ = qbondangle(D23,D24,D34)
			cω,sω = qtorsionangle(D12,D13,D14,D23,D24,D34)
			a = sθ*cω
			b = sθ*sω
			c = -cθ*sω
			d = cθ*cω
			if l==virtualPos
				Q = qprod(Q_before,a,b,c,d)
				lastpos = pos
				break
			else
				Q_virtual = qprod(Q_before,a,b,c,d)
				qmol = rotopt(Q_virtual,D34)
				# TODO: (Emerson) não conseguimos fazer o calculo abaixo da mesma forma que o classicBP?
				vx = qmol.v1 + mol.atoms[virtualLastPos].x-mol.atoms[virtualPos].x
				vy = qmol.v2 + mol.atoms[virtualLastPos].y-mol.atoms[virtualPos].y
				vz = qmol.v3 + mol.atoms[virtualLastPos].z-mol.atoms[virtualPos].z
				if vx*vx + vy*vy + vz*vz > virtual_ε²
					Q_before = qprod(Q_before,a,-b,-c,d)
				else
					Q_before = Q_virtual
				end
				pos = pos+1		
			end
		end
		qmol = rotopt(Q,D34)
		mol.atoms[l].element = NMRdata.info[l,:].nzval[1].atom1
		mol.atoms[l].x = qmol.v1 + mol.atoms[virtualLastPos].x
		mol.atoms[l].y = qmol.v2 + mol.atoms[virtualLastPos].y
		mol.atoms[l].z = qmol.v3 + mol.atoms[virtualLastPos].z
		
		if pruningtest(mol,l,NMRdata,ε)
			if l<n
				quaternionBP_closure(l+1,pos+1,mol,Q)
			else
				nsol=nsol+1
				storage_mol[nsol] = copy(mol)
				return
			end
		end
		if !allmol && nsol>0
			return
		end
		Q = qprod(Q_before,a,-b,-c,d)
		qmol = rotopt(Q,D34)
		mol.atoms[l].x = qmol.v1 + mol.atoms[virtualLastPos].x
		mol.atoms[l].y = qmol.v2 + mol.atoms[virtualLastPos].y
		mol.atoms[l].z = qmol.v3 + mol.atoms[virtualLastPos].z


		if pruningtest(mol,l,NMRdata,ε)  
			if l<n
				quaternionBP_closure(l+1,pos+1,mol,Q)
			else
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				return
			end
		end
	end #closure

	__l, __pos, __mol, __Q = initialization()
	quaternionBP_closure(__l, __pos, __mol, __Q)

	return nsol, storage_mol

end #solver quaternionBP 

