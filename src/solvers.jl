
function quaternion_bp(n :: Int,
			D :: Array{Float64,2},
			nad :: Int,
			ε :: Float64,
			allmol :: Bool)
	#defining bp using closure in order to count
	function qbp(i::Int,n::Int,
		mol::MoleculeType,
		Q::Vector{Quaternion},
		D::Array{Float64,2},
		ε::Float64,
		allmol::Bool)
		if i == 1
			#first atom		
			mol.atoms[1].x = 0.0
			mol.atoms[1].y = 0.0
			mol.atoms[1].z = 0.0
			Q[1] = Quaternion(1.0)
			#second atom
			mol.atoms[2].x = -D[1,2]
			mol.atoms[2].y = 0.0
			mol.atoms[2].z = 0.0
			Q[2] = Quaternion(0.0,0.0,-1.0,0.0)
			# tird atom
			cθ,sθ = qbondangle(3,D)
			Q[3] = Quaternion(sθ,0.0,0.0,cθ)
			Q[3] = Q[2]*Q[3]
			#qmol = Q[3]*Quaternion(0.0,D[3,2],0.0,0.0)*conj(Q[3])
			qmol = rot(Q[3],D[3,2])
			mol.atoms[3].x = qmol.v1 + mol.atoms[2].x
			mol.atoms[3].y = qmol.v2 + mol.atoms[2].y
			mol.atoms[3].z = qmol.v3 + mol.atoms[2].z
			i = 4 # branching starts at atom 4
		end
		λ = 1
		ρ = 1
		cθ,sθ = qbondangle(i,D)
		cω,sω = qtorsionangle(i,D)
		a = sθ*cω
		b = sθ*sω
		c = -cθ*sω
		d = cθ*cω
	#	Q[i] = Quaternion(a,b,c,d)
		Q[i] = Q[i-1]*Quaternion(a,b,c,d)
		qmol = rot(Q[i],D[i,i-1])
		mol.atoms[i].x = qmol.v1 + mol.atoms[i-1].x
		mol.atoms[i].y = qmol.v2 + mol.atoms[i-1].y
		mol.atoms[i].z = qmol.v3 + mol.atoms[i-1].z
		λ  = pruningtest(mol,i,D,ε)
		if λ == 1 
			if i<n
				qbp(i+1,n,mol,Q,D,ε,allmol)
			else
				#vsol[k]=sol

				nsol=nsol+1
				storage_mol[nsol] = copy(mol)
				with_logger(quaternion_bp_logger) do
					@info "Rank n was reached, a solution was found " mol LDE(mol,D,n,nad)
				end
				return 0
			end
		end
		# if only one solution is required, bp stops as soon as the first one is found
		if allmol==false && nsol>0

			#@info "LDE = " LDE(mol,D,n,nad)
			@goto exit
		end
#		@info "partial solution in λ =$(λ) " sol
		
		
		Q[i] = Q[i-1]*Quaternion(a,-b,-c,d)
		qmol = rot(Q[i],D[i,i-1])
		#qmol = Q[i]*Quaternion(0.0,D[i,i-1],0.0,0.0)*conj(Q[i])
		mol.atoms[i].x = qmol.v1 + mol.atoms[i-1].x
		mol.atoms[i].y = qmol.v2 + mol.atoms[i-1].y
		mol.atoms[i].z = qmol.v3 + mol.atoms[i-1].z
		ρ  = pruningtest(mol,i,D,ε)
#		@info "partial solution in ρ =$(ρ) " sol
		if ρ == 1 
			if i<n
				qbp(i+1,n,mol,Q,D,ε,allmol)
			else
				#vsol[k]=sol
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				with_logger(quaternion_bp_logger) do
					@info "Rank n was reached, a solution was found " mol LDE(mol,D,n,nad)
				end
				return 0
			end
		end

	@label exit
	return 0
	end

	mol = MoleculeType(Vector{AtomType}(undef,n),0.0)
	for i=1:n
		mol.atoms[i] = AtomType(0.0,0.0,0.0)
	end
	Q = Vector{Quaternion}(undef,n)
	for i=1:n
			Q[i] = Quaternion(0.0)
	end
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()

	qbp(1,n,mol,Q,D,ε,allmol)

	with_logger(quaternion_bp_logger) do
		@info "number of solutions " nsol
	end
#	for i=1:nsol
#		println(LDE(storage_mol[i],D,n,nad))
		#storage_mol[i].lde = LDE(storage_mol[i],D,n,nad)
#		println(storage_mol[i].atoms[n])
#	end
	return nsol, storage_mol
end



function classical_bp(n :: Int,
			D :: Array{Float64,2},
			nad :: Int,
			ε :: Float64,
			allmol :: Bool)
	#defining bp using closure in order to count
	function bp(i::Int,n::Int,
		mol::MoleculeType,
		C::Vector{Array{Float64,2}},
		D::Array{Float64,2},
		ε::Float64,
		allmol::Bool)
		if i == 1
			#first atom
			mol.atoms[1].x = 0.0
			mol.atoms[1].y = 0.0
			mol.atoms[1].z = 0.0
			C[1] = Diagonal{Float64}(I,4)
			#second atom
			mol.atoms[2].x = -D[1,2]
			mol.atoms[2].y = 0.0
			mol.atoms[2].z = 0.0
			C[2] = zeros(4,4)
			C[2][1,1] = -1.0
			C[2][2,2] = 1.0
			C[2][3,3] = -1.0
			C[2][4,4] = 1.0
			C[2][1,4] = -D[2,1]
			# tird atom
			bdangle = bondangle(3,D)
			mol.atoms[3].x = -D[1,2]+D[2,3]*bdangle.cosθ
			mol.atoms[3].y = D[2,3]*bdangle.sinθ
			mol.atoms[3].z = 0.0
			B = zeros(4,4)
			B[1,1] = -bdangle.cosθ
			B[1,2] = -bdangle.sinθ
			B[1,4] = -D[3,2]*bdangle.cosθ
			B[2,1] = bdangle.sinθ
			B[2,2] = -bdangle.cosθ
			B[2,4] = D[3,2]*bdangle.sinθ
			B[3,3] = 1.0
			B[4,4] = 1.0
			C[3] = C[2]*B
			i = 4 # branching starts at atom 4
		end
		λ = 1
		ρ = 1
		C[i] = prodmatrix(C[i-1],torsionmatrix(i,D,'+'))
		mol.atoms[i].x = C[i][1,4]
		mol.atoms[i].y = C[i][2,4]
		mol.atoms[i].z = C[i][3,4]
		λ  = pruningtest(mol,i,D,ε)
		if λ == 1 
			if i<n
				bp(i+1,n,mol,C,D,ε,allmol)
			else
				#vsol[k]=sol

				nsol=nsol+1
				storage_mol[nsol] = copy(mol)
				with_logger(classical_bp_logger) do
					@info "Rank n was reached, a solution was found " mol LDE(mol,D,n,nad)
				end
				return 0
			end
		end
		# if only one solution is required, bp stops as soon as the first one is found
		if allmol==false && nsol>0

			#@info "LDE = " LDE(mol,D,n,nad)
			@goto exit
		end
#		@info "partial solution in λ =$(λ) " sol
		C[i] = prodmatrix(C[i-1],torsionmatrix(i,D,'-'))
		mol.atoms[i].x = C[i][1,4]
		mol.atoms[i].y = C[i][2,4]
		mol.atoms[i].z = C[i][3,4]
		ρ  = pruningtest(mol,i,D,ε)
#		@info "partial solution in ρ =$(ρ) " sol
		if ρ == 1 
			if i<n
				bp(i+1,n,mol,C,D,ε,allmol)
			else
				#vsol[k]=sol
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				with_logger(classical_bp_logger) do
					@info "Rank n was reached, a solution was found " mol LDE(mol,D,n,nad)
				end
				return 0
			end
		end

	@label exit
	return 0
	end

	mol = MoleculeType(Vector{AtomType}(undef,n),0.0)
	for i=1:n
		mol.atoms[i] = AtomType(0.0,0.0,0.0)
	end
	C = Vector{Array{Float64,2}}(undef,n)
	for i=1:n
		C[i] = zeros(4,4)
	end
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()

	bp(1,n,mol,C,D,ε,allmol)

	with_logger(classical_bp_logger) do
		@info "number of solutions " nsol
	end
#	for i=1:nsol
#		println(LDE(storage_mol[i],D,n,nad))
		#storage_mol[i].lde = LDE(storage_mol[i],D,n,nad)
#		println(storage_mol[i].atoms[n])
#	end
	return nsol, storage_mol
end
