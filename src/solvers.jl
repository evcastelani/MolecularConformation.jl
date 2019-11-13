
function classicBP(NMRdata :: NMRtype,
		   virtual_path :: Vector,
		   ε :: Float64,
		   allmol :: Bool)
	# defining closure
	function classicBP_closure(i :: Int,
				   n :: Int,
				   mol :: MoleculeType,
				   C :: Array{Float64,2})
		if l == 1
			# first atom
			mol.atoms[1].x = 0.0
			mol.atoms[1].y = 0.0
			mol.atoms[1].z = 0.0
			C[1] = Diagonal{Float64}(I,4)
			#second atom
			mol.atoms[2].x = -NMRdata.upperbound[1]
			mol.atoms[2].y = 0.0
			mol.atoms[2].z = 0.0
			C[2] = zeros(4,4)
			C[2][1,1] = -1.0
			C[2][2,2] = 1.0
			C[2][3,3] = -1.0
			C[2][4,4] = 1.0
			C[2][1,4] = -NMRdata.upperbound[1]
			# tird atom
			D12 = NMRdata.upperbound[1]
			D13 = NMRdata.upperbound[2]
			D23 = NMRdata.upperbound[3]
			bdangle = bondangle(D12,D13,D23)
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


	end
end


# esta mantido apenas para cópia de partes do código.
#
function classical_bp(n :: Int,
			D :: Array{Float64,2},
			nad :: Int,
			v:: Vector{Int64},
			optimize::Bool,
			ε :: Float64,
			allmol :: Bool,ndiag::Int)
	#defining bp using closure in order to count
	function bp(i::Int,n::Int,
		mol::MoleculeType,
		C::Vector{Array{Float64,2}},
		D::Array{Float64,2},
		ε::Float64,
		allmol::Bool,ndiag::Int)
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
			repsol = 0
		end
		λ = 1
		ρ = 1

		C[i] = prodmatrix(C[i-1],torsionmatrix(i,D,'+'))
		if optimize==true && v[i]<i-repsol
			println("repeated solution in $(v[i]), level $i " )
			mol.atoms[i].x = mol.atoms[v[i]].x
			mol.atoms[i].y = mol.atoms[v[i]].y
			mol.atoms[i].z = mol.atoms[v[i]].z
			repsol = repsol +1
			λ = 1
		else
			mol.atoms[i].x = C[i][1,4]
			mol.atoms[i].y = C[i][2,4]
			mol.atoms[i].z = C[i][3,4]
			λ  = pruningtest(mol,i,D,ε,ndiag)
		end
		if λ == 1 
			if i<n
				if last_level < i
					last_level = i
					last_level_count = 1
				end
				if last_level == i
					last_level_count = last_level_count +1
					println("level $(last_level), count $(last_level_count)")
				end
#				println("level $i, sol = [$(mol.atoms[i].x),$(mol.atoms[i].y),$(mol.atoms[i].z)]")
				bp(i+1,n,mol,C,D,ε,allmol,ndiag)
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
		if optimize==true && v[i]<i-repsol 
		
			repsol = repsol+1
			ρ = 1
			println("repeated solution  $(v[i]), level $i " )
			mol.atoms[i].x = mol.atoms[v[i]].x
			mol.atoms[i].y = mol.atoms[v[i]].y
			mol.atoms[i].z = mol.atoms[v[i]].z
		else
			mol.atoms[i].x = C[i][1,4]
			mol.atoms[i].y = C[i][2,4]
			mol.atoms[i].z = C[i][3,4]
			ρ  = pruningtest(mol,i,D,ε,ndiag)
		end
#		@info "partial solution in ρ =$(ρ) " sol
		if ρ == 1 
			if i<n
				if last_level < i
					last_level = i
					last_level_count = 1
				end
				if last_level == i
					last_level_count = last_level_count +1
					
					println("level $(last_level), count $(last_level_count)")
				end

#				println("level $i, sol = [$(mol.atoms[i].x),$(mol.atoms[i].y),$(mol.atoms[i].z)]")
				bp(i+1,n,mol,C,D,ε,allmol,ndiag)
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
	repsol = 0
	storage_mol = Dict{Int64,MoleculeType}()
	last_level = 1
	last_level_count = 1
	bp(1,n,mol,C,D,ε,allmol,ndiag)
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
