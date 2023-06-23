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
		allmol :: Bool, 
		time_limit::Union{Period,Nothing}=nothing;
		output_to_symBP::Bool=false)
	
	if time_limit !== nothing
		start =  Dates.now()
		time_elapsed = Second(0.0)
	end

	n = NMRdata.dim
	if n < 3
		ArgumentError("Invalid dimension of NMRdata")
	end

	virtual_ε² = virtual_ε*virtual_ε 
	ϵₛ = ε
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()

	if !allmol && output_to_symBP
		path_to_sol = zeros(Bool,n)
	end

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
		mol.atoms[3].x = D23*cθ-D12
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

		if time_limit !== nothing
			time_elapsed = Dates.now()-start
			if time_elapsed>time_limit && l<n
				error("Time limit reached without found a solution!")
			end
		end

		C_before = copy(C)
		B = Array{Float64,2}(undef,4,4)
		while true
			virtualPos = NMRdata.virtual_path[pos]
			virtualLastPos = NMRdata.virtual_path[pos-1]
			
			D14 = NMRdata.info[NMRdata.virtual_path[pos-3],virtualPos].dist
			D24 = NMRdata.info[NMRdata.virtual_path[pos-2],virtualPos].dist
			D34 = NMRdata.info[virtualLastPos,virtualPos].dist
			D12 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos-2]].dist
			D13 = NMRdata.info[NMRdata.virtual_path[pos-3],virtualLastPos].dist
			D23 = NMRdata.info[NMRdata.virtual_path[pos-2],virtualLastPos].dist
			
			cθ,sθ = bondangle(D23,D24,D34)
			cω,sω = torsionangle(D12,D13,D14,D23,D24,D34)
			torsionmatrix(B,cθ,sθ,cω,sω,D34)
			if l==virtualPos
				prodmatrix(C,C_before,B)
				break
			else
				cpx = mol.atoms[virtualPos].x
				cpy = mol.atoms[virtualPos].y
				cpz = mol.atoms[virtualPos].z

				Virtual_Torsion = prodmatrix(C_before,B)

				if (Virtual_Torsion[1,4]- cpx)^2+(Virtual_Torsion[2,4]- cpy)^2+(Virtual_Torsion[3,4]- cpz)^2> virtual_ε²
					reflectmatrix(B)
					C_before = prodmatrix(C_before,B)
				else
					C_before = Virtual_Torsion
				end
				pos = pos+1		
			end
		end

		#mol.atoms[l].element = NMRdata.info[l,:].nzval[1].atom1 we not need anymore
		changeposition(mol.atoms[l], C[1,4], C[2,4], C[3,4])
		
		if pruningtest(mol,l,NMRdata,ϵₛ) 
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
		
		reflectmatrix(B)
		prodmatrix(C,C_before,B)
		changeposition(mol.atoms[l], C[1,4], C[2,4], C[3,4])
		
		if pruningtest(mol,l,NMRdata,ϵₛ) #preciso modificar 
			if !allmol && output_to_symBP
				path_to_sol[l] = true
			end
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
	
	if !allmol && output_to_symBP
		return ConformationOutput(classicBP , nsol, storage_mol), BitVector(path_to_sol)
	end
	return ConformationOutput(classicBP , nsol, storage_mol)

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
		allmol :: Bool, 
		time_limit::Union{Period,Nothing}=nothing;
		output_to_symBP::Bool=false)
	
	if time_limit !== nothing
		start =  Dates.now()
		time_elapsed = Second(0.0)
	end

	n = NMRdata.dim
	if n < 3
		ArgumentError("Invalid dimension of NMRdata")
	end

	virtual_ε² = virtual_ε*virtual_ε 
	ϵₛ = ε
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()
	
	if !allmol && output_to_symBP
		path_to_sol = zeros(Bool,n)
	end

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
		# cθ,sθ = bondangle(D12,D13,D23)
		# cθ = sqrt(0.5 + 0.5*cθ)
		# sθ = sqrt(0.5 - 0.5*sθ)
		cθ,sθ = qbondangle(D12,D13,D23)
		Q = Quaternion(0.0,-cθ,-sθ,0.0)
		d = 2.0*D23
		qmol = Quaternion(0.0,d*(cθ*cθ-0.5),d*(cθ*sθ),0.0)
		mol.atoms[3].element = NMRdata.info[3,:].nzval[1].atom1
		mol.atoms[3].x = qmol.v1 + mol.atoms[2].x
		mol.atoms[3].y = qmol.v2 + mol.atoms[2].y
		mol.atoms[3].z = qmol.v3 + mol.atoms[2].z
		# cost: 4 + 5c_m
		# qcos(θ) cost = 3 + 5c_m + c_d + c_√
		# qsin(θ) cost = 1 + c_√
		# cost total: 8 + 10c_m + c_d + 2c_√
		return 4,4,mol,Q
	end
	
	# defining closure
	function quaternionBP_closure(l :: Int64,
									pos::Int64,
									mol :: MoleculeType,
									Q :: Quaternion)

		if time_limit !== nothing
			time_elapsed = Dates.now()-start
			if time_elapsed>time_limit && l<n
				error("Time limit reached without found a solution!")
			end
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

			D14 = NMRdata.info[NMRdata.virtual_path[pos-3],virtualPos].dist
			D24 = NMRdata.info[NMRdata.virtual_path[pos-2],virtualPos].dist
			D34 = NMRdata.info[virtualLastPos,virtualPos].dist
			D12 = NMRdata.info[NMRdata.virtual_path[pos-3],NMRdata.virtual_path[pos-2]].dist
			D13 = NMRdata.info[NMRdata.virtual_path[pos-3],virtualLastPos].dist
			D23 = NMRdata.info[NMRdata.virtual_path[pos-2],virtualLastPos].dist
			
			cθ,sθ = qbondangle(D23,D24,D34)
			cω,sω = qtorsionangle(D12,D13,D14,D23,D24,D34)
			a = sθ*cω
			b = sθ*sω
			c = -cθ*sω
			d = cθ*cω
			# q_i^0 cost: 4c_m
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
		#mol.atoms[l].element = NMRdata.info[l,:].nzval[1].atom1 we not need anymore
		changeposition(mol.atoms[l], qmol.v1 + mol.atoms[virtualLastPos].x, 
									 qmol.v2 + mol.atoms[virtualLastPos].y,
									 qmol.v3 + mol.atoms[virtualLastPos].z)

		if pruningtest(mol,l,NMRdata,ϵₛ)
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

		if pruningtest(mol,l,NMRdata,ϵₛ)  
			if !allmol && output_to_symBP
				path_to_sol[l] = true
			end
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

	if !allmol && output_to_symBP
		return ConformationOutput(quaternionBP, nsol, storage_mol), BitVector(path_to_sol)
	end
	return ConformationOutput(quaternionBP, nsol, storage_mol)
end #solver quaternionBP 


