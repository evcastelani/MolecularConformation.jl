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
	function classicBP_closure(l :: Int,
				   pos::Int,
				   mol :: MoleculeType,
				   C :: Array{Float64,2})
		if l == 1
			# first atom
			mol.atoms[1].x = 0.0
			mol.atoms[1].y = 0.0
			mol.atoms[1].z = 0.0
			#second atom
			mol.atoms[2].x = -NMRdata.info[1,2].dist
			mol.atoms[2].y = 0.0
			mol.atoms[2].z = 0.0
			C = zeros(4,4)
			C[1,1] = -1.0
			C[2,2] = 1.0
			C[3,3] = -1.0
			C[4,4] = 1.0
			C[1,4] = -NMRdata.info[1,2].dist
			# tird atom
			D12 = NMRdata.info[1,2].dist
			D13 = NMRdata.info[1,3].dist
			D23 = NMRdata.info[2,3].dist
			D14 = 0.0
			D24 = 0.0
			D34 = 0.0
			cθ,sθ = bondangle(D12,D13,D23)
			cω,sω = (0.0,0.0)
			mol.atoms[3].x = -D12+D23*cθ
			mol.atoms[3].y = D23*sθ
			mol.atoms[3].z = 0.0
			B = zeros(4,4)
			B[1,1] = -cθ
			B[1,2] = -sθ
			B[1,4] = -D23*cθ
			B[2,1] = sθ
			B[2,2] = -cθ
			B[2,4] = D23*sθ
			B[3,3] = 1.0
			B[4,4] = 1.0
			C = C*B
			l = 4 # branching starts at atom 4
			pos = 4 # position in virtual path
		end

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
			cω,sω = torsionangle(D12,D13,D14,D23,D24,D34)
			@debug "l value = $(l) and NMRdatavalue = $(NMRdata.virtual_path[pos]) in position $(pos)"
			if l==NMRdata.virtual_path[pos]
				B = torsionmatrix(cθ,sθ,cω,sω,D34,true)
				C = C_before*B
				keep = false
			else
				B = torsionmatrix(cθ,sθ,cω,sω,D34,true)
				cpx = mol.atoms[NMRdata.virtual_path[pos]].x
				cpy = mol.atoms[NMRdata.virtual_path[pos]].y
				cpz = mol.atoms[NMRdata.virtual_path[pos]].z

				Virtual_Torsion = C_before*B
				if sqrt((Virtual_Torsion[1,4]- cpx)^2+(Virtual_Torsion[2,4]- cpy)^2+(Virtual_Torsion[3,4]- cpz)^2)> virtual_ε
					B = torsionmatrix(cθ,sθ,cω,sω,D34,false)
				end
				C_before = C_before*B	


				@debug "virtual atom position  " C_before[1,4],C_before[2,4],C_before[3,4]
				pos = pos+1		
			end
		end
				
		mol.atoms[l].x = C[1,4]
		mol.atoms[l].y = C[2,4]
		mol.atoms[l].z = C[3,4]
		λ  = pruningtest(mol,l,NMRdata,ε) 
		@debug "C at level $(l) right side " C
		if λ == 1 
			if l<n
				@debug "Partial solution by right side at level $(l)"  mol
				classicBP_closure(l+1,pos+1,mol,C)
			else
				nsol=nsol+1
				storage_mol[nsol] = copy(mol)
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
		end
		if allmol==false && nsol>0
			@debug "number of solutions"  nsol
			#@info "LDE = " LDE(mol,D,n,nad)
			@goto exit
		end
		B = torsionmatrix(cθ,sθ,cω,sω,D34,false)
		C = C_before*B
		mol.atoms[l].x = C[1,4]
		mol.atoms[l].y = C[2,4]
		mol.atoms[l].z = C[3,4]
		ρ  = pruningtest(mol,l,NMRdata,ε) #preciso modificar
		@debug "C at level $(l) left side " C
		if ρ == 1 
			if l<n
				@debug "Partial solution by left side at level $(l)" mol
				classicBP_closure(l+1,pos+1,mol,C)
			else
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
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
	classicBP_closure(1,1,mol,C)
	return nsol, storage_mol

end #solver classicBP

"""
```
quaternionBP :: Function
```
This function defines a new solver. The implementation follows ideas describing in:

Fidalgo, F. Using Quaternion Geometric Algebra for efficient rotiations in Branch and Prune Algorithtm to solve the Discretizable Molecular Distance Geometry Problem. In: Proceedings of AGACSE 2018, Campinas-SP, Brazil.
"""
function quaternionBP(NMRdata :: NMRType,
		   ε :: Float64,
		   virtual_ε :: Float64,
		   allmol :: Bool)
	# defining closure
	function quaternionBP_closure(l :: Int,
				   pos::Int,
				   mol :: MoleculeType,
				   Q :: Quaternion)
		if l == 1
			# first atom
			mol.atoms[1].x = 0.0
			mol.atoms[1].y = 0.0
			mol.atoms[1].z = 0.0
			#second atom
			mol.atoms[2].x = -NMRdata.info[1,2].dist
			mol.atoms[2].y = 0.0
			mol.atoms[2].z = 0.0
			Q_before = Quaternion(0.0,0.0,-1.0,0.0)
			# tird atom
			D12 = NMRdata.info[1,2].dist
			D13 = NMRdata.info[1,3].dist
			D23 = NMRdata.info[2,3].dist
			D14 = 0.0
			D24 = 0.0
			D34 = 0.0
			cθ,sθ = qbondangle(D12,D13,D23)
			Q = qprod(Q_before,Quaternion(sθ,0.0,0.0,cθ))
			qmol = rot(Q,D23)
			mol.atoms[3].x = qmol.v1 + mol.atoms[2].x
			mol.atoms[3].y = qmol.v2 + mol.atoms[2].y
			mol.atoms[3].z = qmol.v3 + mol.atoms[2].z
			l = 4 # branching starts at atom 4
			pos = 4 # position in virtual path
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
			cω,sω = qtorsionangle(D12,D13,D14,D23,D24,D34)
			a = sθ*cω
			b = sθ*sω
			c = -cθ*sω
			d = cθ*cω
			@debug "l value = $(l) and NMRdatavalue = $(NMRdata.virtual_path[pos]) in position $(pos)"
			if l==NMRdata.virtual_path[pos]
				Q = qprod(Q_before,Quaternion(a,b,c,d))
				keep = false
				lastpos = pos
			else
				Q_virtual = qprod(Q_before,Quaternion(a,b,c,d))
				qmol = rot(Q_virtual,D34) #? is 34?
				vx = qmol.v1 + mol.atoms[NMRdata.virtual_path[pos-1]].x-mol.atoms[NMRdata.virtual_path[pos]].x
				vy = qmol.v2 + mol.atoms[NMRdata.virtual_path[pos-1]].y-mol.atoms[NMRdata.virtual_path[pos]].y
				vz = qmol.v3 + mol.atoms[NMRdata.virtual_path[pos-1]].z-mol.atoms[NMRdata.virtual_path[pos]].z

				if sqrt(vx^2+vy^2+vz^2)> virtual_ε
					Q_before = qprod(Q_before,Quaternion(a,-b,-c,d))
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
		qmol = rot(Q,D34)
		mol.atoms[l].x = qmol.v1 + mol.atoms[NMRdata.virtual_path[lastpos-1]].x
		mol.atoms[l].y = qmol.v2 + mol.atoms[NMRdata.virtual_path[lastpos-1]].y
		mol.atoms[l].z = qmol.v3 + mol.atoms[NMRdata.virtual_path[lastpos-1]].z
		@debug "candidate atom by right side at level $(l) = $(mol.atoms[l].x),$(mol.atoms[l].y),$(mol.atoms[l].z)"
		λ  = pruningtest(mol,l,NMRdata,ε) 
		if λ == 1 
			if l<n
				@debug "Partial solution by right side at level $(l)"  mol
				quaternionBP_closure(l+1,pos+1,mol,Q)
			else
				nsol=nsol+1
				storage_mol[nsol] = copy(mol)
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
		end
		if allmol==false && nsol>0
			@debug "number of solutions"  nsol
			@goto exit
		end
		Q = qprod(Q_before,Quaternion(a,-b,-c,d))
		qmol = rot(Q,D34)
		mol.atoms[l].x = qmol.v1 + mol.atoms[NMRdata.virtual_path[lastpos-1]].x
		mol.atoms[l].y = qmol.v2 + mol.atoms[NMRdata.virtual_path[lastpos-1]].y
		mol.atoms[l].z = qmol.v3 + mol.atoms[NMRdata.virtual_path[lastpos-1]].z
		@debug "candidate atom by left side at level $(l) = $(mol.atoms[l].x),$(mol.atoms[l].y),$(mol.atoms[l].z)"
		ρ  = pruningtest(mol,l,NMRdata,ε) #preciso modificar
		if ρ == 1 
			if l<n
				@debug "Partial solution by left side at level $(l)" mol
				quaternionBP_closure(l+1,pos+1,mol,Q)
			else
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
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
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()
	quaternionBP_closure(1,1,mol,Q)

	return nsol, storage_mol

end #solver quaternionBO 

