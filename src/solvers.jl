
function classicBP(NMRdata :: NMRType,
		   ε :: Float64,
		   allmol :: Bool)
	# defining closure
	function classicBP_closure(l :: Int,
				   pos::Int,
				   mol :: MoleculeType,
				   sign:: Vector{Bool}, 
				   C :: Vector{Array{Float64,2}})
		if l == 1
			# first atom
			mol.atoms[1].x = 0.0
			mol.atoms[1].y = 0.0
			mol.atoms[1].z = 0.0
			sign[1] = true
			C[1] = Diagonal{Float64}(I,4)
			#second atom
			#mol.atoms[2].x = -D[1,2]
			mol.atoms[2].x = -NMRdata.info[1,2].dist
			mol.atoms[2].y = 0.0
			mol.atoms[2].z = 0.0
			sign[2] = true
			C[2] = zeros(4,4)
			C[2][1,1] = -1.0
			C[2][2,2] = 1.0
			C[2][3,3] = -1.0
			C[2][4,4] = 1.0
			C[2][1,4] = -NMRdata.info[1,2].dist
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
			sign[3] = true
			B = zeros(4,4)
			B[1,1] = -cθ
			B[1,2] = -sθ
			B[1,4] = -D23*cθ
			B[2,1] = sθ
			B[2,2] = -cθ
			B[2,4] = D23*sθ
			B[3,3] = 1.0
			B[4,4] = 1.0
			C[3] = C[2]*B
			l = 4 # branching starts at atom 4
			pos = 4 # position in virtual path
		end

		λ = 1
		ρ = 1
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
				C[l] = C[l-1]*B
				sign[l] = true
				keep = false
			else
				B = torsionmatrix(cθ,sθ,cω,sω,D34,true)
				Aux = zeros(4,4)
				copyto!(Aux,C[l-1]*B)
				cpx = mol.atoms[NMRdata.virtual_path[pos]].x
				cpy = mol.atoms[NMRdata.virtual_path[pos]].y
				cpz = mol.atoms[NMRdata.virtual_path[pos]].z
				mol.atoms[NMRdata.virtual_path[pos]].x = Aux[1,4]
				mol.atoms[NMRdata.virtual_path[pos]].y = Aux[2,4]
				mol.atoms[NMRdata.virtual_path[pos]].z = Aux[3,4]
				λ  = pruningtest(mol,NMRdata.virtual_path[pos],NMRdata,ε)
				if 	λ == 0
					B = torsionmatrix(cθ,sθ,cω,sω,D34,false)
					Aux = C[l-1]*B
				end
				C[l-1] = Aux

				mol.atoms[NMRdata.virtual_path[pos]].x = cpx
				mol.atoms[NMRdata.virtual_path[pos]].y = cpy
				mol.atoms[NMRdata.virtual_path[pos]].z = cpz
				pos = pos+1		
			end
		end
				
		mol.atoms[l].x = C[l][1,4]
		mol.atoms[l].y = C[l][2,4]
		mol.atoms[l].z = C[l][3,4]
		λ  = pruningtest(mol,l,NMRdata,ε) 
		if λ == 1 
			if l<n
				@debug "Partial solution by right side at level $(l)" mol
				classicBP_closure(l+1,pos+1,mol,sign,C)
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
		C[l] = C[l-1]*B
		mol.atoms[l].x = C[l][1,4]
		mol.atoms[l].y = C[l][2,4]
		mol.atoms[l].z = C[l][3,4]
		sign[l]=false
		ρ  = pruningtest(mol,l,NMRdata,ε) #preciso modificar
		if ρ == 1 
			if l<n
				@debug "Partial solution by left side at level $(l) " mol
				classicBP_closure(l+1,pos+1,mol,sign,C)
			else
				nsol = nsol+1
				storage_mol[nsol] = copy(mol)				
				@debug "Rank n was reached, a solution was found " 
				return 0
			end
		end

		@label exit
		return 0
	end
	# end of bp_closure
	n = NMRdata.dim
	mol = MoleculeType(Vector{AtomType}(undef,n),0.0)
	for i=1:n
		mol.atoms[i] = AtomType(0.0,0.0,0.0)
	end
	C = Vector{Array{Float64,2}}(undef,n)
	signal = Vector{Bool}(undef,n)
	for i=1:n
		signal[i] = true
		C[i] = zeros(4,4)
	end
	nsol = 0
	storage_mol = Dict{Int64,MoleculeType}()
	classicBP_closure(1,1,mol,signal,C)
	return nsol, storage_mol

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
