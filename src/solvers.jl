function classicBPseq(NMRdata::NMRType,
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
	C_list[3] = C
	B = Array{Array{Float64,2}}(undef,n) # to acess level l it is need to put l-3
	B[3] = zeros(4,4)
	while l<=n
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
			@debug "l value = $(l) and NMRdatavalue = $(NMRdata.virtual_path[pos]) in position $(pos)"
			B[l] = torsionmatrix(cθ,sθ,cω,sω,D34)
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
				nop_vpath += [24,33,0,0]
				nop_node += [23,33,0,0]

				if sqrt((Virtual_Torsion[1,4]- cpx)^2+(Virtual_Torsion[2,4]- cpy)^2+(Virtual_Torsion[3,4]- cpz)^2)> virtual_ε
					B[l] = torsionmatrix(B[l])
					C_list[l-1] = prodmatrix(C_list[l-1],B[l])
					nop_vpath += [24,33,0,0]
					nop_node += [24,33,0,0]

				else
					#copyto!(C_before,Virtual_Torsion) 
					C_list[l-1] = Virtual_Torsion
				end
				#@debug "virtual atom position  " C_before[1,4],C_before[2,4],C_before[3,4]
				pos = pos+1		
			end
		end


		if explore_right_side[l] == false
			mol.atoms[l].element = NMRdata.info[l,:].nzval[1].atom1		
			mol.atoms[l].x = C_list[l][1,4]
			mol.atoms[l].y = C_list[l][2,4]
			mol.atoms[l].z = C_list[l][3,4]
			count = [0,0,0,0]
			λ , count  = pruningtest(mol,l,NMRdata,ε,count) 
			nop_ddf += count
			@debug "C at level $(l) right side " C[l]
			if λ == 1 
				if l<n
					@debug "Partial solution by right side at level $(l)"  mol
					n_branch +=1
					#classicBP_closure(l+1,pos+1,mol,C)
					l += 1
					explore_right_side[l]=false
				else
					nsol=nsol+1
					storage_mol[nsol] = copy(mol)
					@debug "Rank n was reached, a solution was found " 
					return 0
				end
			else
				n_prune += 1
				explore_right_side[l] = true
			end
		end
		if explore_right_side[l] == true 
			B[l] = torsionmatrix(B[l])
			#nop_node += [0,0,0,0]
			C_list[l] = prodmatrix(C_list[l-1],B[l])# tenho que otimizar este calculo
			nop_node += [24,33,0,0]
			mol.atoms[l].x = C_list[l][1,4]
			mol.atoms[l].y = C_list[l][2,4]
			mol.atoms[l].z = C_list[l][3,4]
			count = [0,0,0,0]
			ρ ,count = pruningtest(mol,l,NMRdata,ε,count) #preciso modificar
			nop_ddf += count 
			@debug "C at level $(l) left side " C
			if ρ == 1 
				if l<n
					@debug "Partial solution by left side at level $(l)" mol
					n_branch += 1
					l += 1
					explore_right_side[l]=false
					#classicBP_closure(l+1,pos+1,mol,C)
				else
					nsol = nsol+1
					storage_mol[nsol] = copy(mol)				
					@debug "Rank n was reached, a solution was found " 
					explore_right_side = zeros(Bool,n)
				end
			else
				k = l-1
				while explore_right_side[k] == true
					k -= 1
				end
				explore_right_side[k] = true
				n_prune += 1
			end


		end

	end
	return nsol, storage_mol,Counter(nop_node,nop_vpath,nop_ddf,n_branch,n_prune)
end
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
		allmol :: Bool)
	# defining closure
	function classicBP_closure(l :: Int64,
			pos::Int64,
			mol :: MoleculeType,
			C :: Array{Float64,2})
		if l == 1
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
		end
		B = zeros(4,4)
		λ = 1
		ρ = 1
		C_before = zeros(4,4)
		copyto!(C_before,C)
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
			@debug "l value = $(l) and NMRdatavalue = $(NMRdata.virtual_path[pos]) in position $(pos)"
			B = torsionmatrix(cθ,sθ,cω,sω,D34,B,true)
			if l==NMRdata.virtual_path[pos]
				nop_node += [0,7,0,0] #torsion matrix
				nop_node += [3,6,1,1] # bond angle
				nop_node += [10,20,4,2] # bad torsion angle
				C = prodmatrix(C_before,B)
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

				Virtual_Torsion = prodmatrix(C_before,B)
				nop_vpath += [24,33,0,0]
				nop_node += [23,33,0,0]

				if sqrt((Virtual_Torsion[1,4]- cpx)^2+(Virtual_Torsion[2,4]- cpy)^2+(Virtual_Torsion[3,4]- cpz)^2)> virtual_ε
					B = torsionmatrix(cθ,sθ,cω,sω,D34,B,false)
					C_before = prodmatrix(C_before,B)
					nop_vpath += [24,33,0,0]
					nop_node += [24,33,0,0]

				else
					copyto!(C_before,Virtual_Torsion) 
				end
				@debug "virtual atom position  " C_before[1,4],C_before[2,4],C_before[3,4]
				pos = pos+1		
			end
		end

		mol.atoms[l].element = NMRdata.info[l,:].nzval[1].atom1		
		mol.atoms[l].x = C[1,4]
		mol.atoms[l].y = C[2,4]
		mol.atoms[l].z = C[3,4]
		count = [0,0,0,0]
		λ , count  = pruningtest(mol,l,NMRdata,ε,count) 
		nop_ddf += count
		@debug "C at level $(l) right side " C
		if λ == 1 
			if l<n
				@debug "Partial solution by right side at level $(l)"  mol
				n_branch +=1
				classicBP_closure(l+1,pos+1,mol,C)
			else
				nsol=nsol+1
				storage_mol[nsol] = copy(mol)
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
		else
			n_prune += 1
		end
		if allmol==false && nsol>0
			@debug "number of solutions"  nsol
			#@info "LDE = " LDE(mol,D,n,nad)
			@goto exit
		end
		B = torsionmatrix(cθ,sθ,cω,sω,D34,B,false)
		#nop_node += [0,0,0,0]
		C = prodmatrix(C_before,B)# tenho que otimizar este calculo
		nop_node += [24,33,0,0]
		mol.atoms[l].x = C[1,4]
		mol.atoms[l].y = C[2,4]
		mol.atoms[l].z = C[3,4]
		count = [0,0,0,0]
		ρ ,count = pruningtest(mol,l,NMRdata,ε,count) #preciso modificar
		nop_ddf += count 
		@debug "C at level $(l) left side " C
		if ρ == 1 
			if l<n
				@debug "Partial solution by left side at level $(l)" mol
				n_branch += 1
				classicBP_closure(l+1,pos+1,mol,C)
			else
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
		else
			n_prune += 1
		end

		@label exit
		return 0
	end # closure
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
	classicBP_closure(1,1,mol,C)
	return nsol, storage_mol,Counter(nop_node,nop_vpath,nop_ddf,n_branch,n_prune)

end #solver classicBP

""" 
```
quaternionBP :: Function
```
This function defines a new solver. The implementation follows ideas describing in:

Fidalgo, F. Using Quaternion Geometric Algebra for efficient rotations in Branch and Prune Algorithtm to solve the Discretizable Molecular Distance Geometry Problem. In: Proceedings of AGACSE 2018, Campinas-SP, Brazil.
"""
#qbondangle = [4,5,2,2]
#qtorsionangle = [11,19]
function quaternionBP(NMRdata :: NMRType,
		ε :: Float64,
		virtual_ε :: Float64,
		allmol :: Bool)
	# defining closure
	function quaternionBP_closure(l :: Int64,
			pos::Int64,
			mol :: MoleculeType,
			Q :: Quaternion)
		if l == 1
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
			#Q_before = Quaternion(0.0,0.0,-1.0,0.0)
			# tird atom
			D12 = NMRdata.info[1,2].dist
			D13 = NMRdata.info[1,3].dist
			D23 = NMRdata.info[2,3].dist
			D14 = 0.0
			D24 = 0.0
			D34 = 0.0
			cθ,sθ = qbondangle(D12,D13,D23)
			nop_node += [4,5,1,2]
			#Q = qprod(Q_before,Quaternion(sθ,0.0,0.0,cθ)) 
			#nop_node += [4,8,0,0]
			#qmol = rotopt(Q,D23)
			#nop_node += [4,10,0,0]
			Q = Quaternion(0.0,-cθ,-sθ,0.0)
			d = 2.0*D23
			qmol = Quaternion(0.0,d*(cθ^2-0.5),d*(cθ*sθ),0.0)
			nop_node = [1,4,0,0]
			mol.atoms[3].element = NMRdata.info[3,:].nzval[1].atom1
			mol.atoms[3].x = qmol.v1 + mol.atoms[2].x
			mol.atoms[3].y = qmol.v2 + mol.atoms[2].y
			mol.atoms[3].z = qmol.v3 + mol.atoms[2].z
			nop_node += [3,0,0,0]
			l = Int64(4) # branching starts at atom 4
			pos = Int64(4) # position in virtual path
		end
		lastpos = 1
		λ = 1
		ρ = 1
		a = 0.0
		b = 0.0
		c = 0.0
		d = 0.0 
		Q_before = Q
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
			cθ,sθ = qbondangle(D23,D24,D34)
			nop_node += [4,5,1,2]
			cω,sω = qtorsionangle(D12,D13,D14,D23,D24,D34)
			nop_node += [11,16,1,3]
			a = sθ*cω
			b = sθ*sω
			c = -cθ*sω
			d = cθ*cω
			nop_node += [0,4,0,0]
			@debug "l value = $(l) and NMRdatavalue = $(NMRdata.virtual_path[pos]) in position $(pos)"
			if l==NMRdata.virtual_path[pos]
				Q = qprod(Q_before,Quaternion(a,b,c,d))
				keep = false
				lastpos = pos
				nop_node += [12,16,0,0]
			else
				nop_vpath += [19,21,2,5] #bond+torsion+product
				Q_virtual = qprod(Q_before,Quaternion(a,b,c,d))
				nop_vpath += [12,16,0,0]
				nop_node += [12,16,0,0]
				qmol = rotopt(Q_virtual,D34)
				nop_vpath += [4,10,0,0]
				nop_node += [4,10,0,0]
				vx = qmol.v1 + mol.atoms[NMRdata.virtual_path[pos-1]].x-mol.atoms[NMRdata.virtual_path[pos]].x
				vy = qmol.v2 + mol.atoms[NMRdata.virtual_path[pos-1]].y-mol.atoms[NMRdata.virtual_path[pos]].y
				vz = qmol.v3 + mol.atoms[NMRdata.virtual_path[pos-1]].z-mol.atoms[NMRdata.virtual_path[pos]].z
				nop_node += [3,0,0,0]
				nop_vpath += [3,0,0,0]
				if sqrt(vx^2+vy^2+vz^2)> virtual_ε
					Q_before = qprod(Q_before,Quaternion(a,-b,-c,d))
					nop_node += [12,16,0,0]
					nop_vpath +=[12,16,0,0]
				else
					Q_before = Q_virtual
				end
				#qmol_aux = rot(Q_before,D34) #? is 34?
				#vx = qmol_aux.v1 + mol.atoms[NMRdata.virtual_path[pos-1]].x
				#vy = qmol_aux.v2 + mol.atoms[NMRdata.virtual_path[pos-1]].y
				#vz = qmol_aux.v3 + mol.atoms[NMRdata.virtual_path[pos-1]].z
				#@debug "repeated atom = $(vx),$(vy),$(vz)"
				pos = pos+1		
			end
		end
		qmol = rotopt(Q,D34)
		nop_node += [4,10,0,0] 
		mol.atoms[l].element = NMRdata.info[l,:].nzval[1].atom1
		mol.atoms[l].x = qmol.v1 + mol.atoms[NMRdata.virtual_path[lastpos-1]].x
		mol.atoms[l].y = qmol.v2 + mol.atoms[NMRdata.virtual_path[lastpos-1]].y
		mol.atoms[l].z = qmol.v3 + mol.atoms[NMRdata.virtual_path[lastpos-1]].z
		nop_node += [3,0,0,0]
		@debug "candidate atom by right side at level $(l) = $(mol.atoms[l].x),$(mol.atoms[l].y),$(mol.atoms[l].z)"
		count = [0,0,0,0]
		λ , count  = pruningtest(mol,l,NMRdata,ε,count) 
		nop_ddf += count 
		if λ == 1 
			if l<n
				@debug "Partial solution by right side at level $(l)"  mol
				n_branch +=1
				quaternionBP_closure(l+1,pos+1,mol,Q)
			else
				nsol=nsol+1
				storage_mol[nsol] = copy(mol)
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
		else
			n_prune += 1
			# println("count prune = $(count_prun)")
		end
		if allmol==false && nsol>0
			@debug "number of solutions"  nsol
			@goto exit
		end
		Q = qprod(Q_before,Quaternion(a,-b,-c,d))
		nop_node += [12,16,0,0] 
		qmol = rotopt(Q,D34)
		nop_node += [4,10,0,0]
		mol.atoms[l].x = qmol.v1 + mol.atoms[NMRdata.virtual_path[lastpos-1]].x
		mol.atoms[l].y = qmol.v2 + mol.atoms[NMRdata.virtual_path[lastpos-1]].y
		mol.atoms[l].z = qmol.v3 + mol.atoms[NMRdata.virtual_path[lastpos-1]].z
		nop_node += [3,0,0,0]
		@debug "candidate atom by left side at level $(l) = $(mol.atoms[l].x),$(mol.atoms[l].y),$(mol.atoms[l].z)"
		count = [0,0,0,0]
		ρ,count  = pruningtest(mol,l,NMRdata,ε,count) 
		nop_ddf += count
		if ρ == 1 
			if l<n
				@debug "Partial solution by left side at level $(l)" mol
				n_branch += 1
				quaternionBP_closure(l+1,pos+1,mol,Q)
			else
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
		else
			n_prune +=1
			# println("count prune = $(count_prun)")
		end

		@label exit
		return 0

	end #closure

	n = NMRdata.dim
	mol = MoleculeType(Vector{AtomType}(undef,n),0.0)
	for i=1:n
		mol.atoms[i] = AtomType(0.0,0.0,0.0)
	end
	Q = Quaternion(0.0,0.0,0.0,0.0)
	nsol = Int64(0)
	storage_mol = Dict{Int64,MoleculeType}()
	n_prune = 0
	n_branch = 0
	#count_nop = [+-,*,/,√]
	nop_node = [0,0,0,0]
	nop_ddf = [0,0,0,0]
	nop_vpath = [0,0,0,0]
	quaternionBP_closure(1,1,mol,Q)
	#println(" *** number of main operations evaluated $(count_op)")
	return nsol, storage_mol,Counter(nop_node,nop_vpath,nop_ddf,n_branch,n_prune)

end #solver quaternionBP 

