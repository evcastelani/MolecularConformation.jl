# these types are mandatory indepedent of used solver##############################
#
"""
```
NMRinfo :: Type
```
This type contains the information for `NMRtype`. In this point of the project it contains two fields: dist, atom1 and atom2. Essentially, dist is the distance between atom1 and atom2 which are strings related to the name of atoms.  
"""
struct NMRInfo
	dist::Float64
	atom1::String
	atom2::String
end

"""
```
NMRType :: Type
```
This type is used as input data. Essentially, this type carry itself fields for internal solvers. Although it is possible to input data directly into the type, it is highly recommended that this type be constructed from the function `preprocessing`, like in the following example 

## Example

```julia-repl
nmr_data_example = preprocessing("1A57.nmr")
```
As return an element of `NMRtype` is created. We are assuming that a `1A57.nmr` file is provided.
"""
struct NMRType
	virtual_path :: Vector{Int64}
	additional_distance::Vector{Vector{Int64}}
	info :: SparseMatrixCSC{NMRInfo,Int64}
	dim :: Int64
end

"""
```
preprocessing :: Function
```
This function is able to read a `.nmr` file and generate a `NMRType`, which is an internal type with fields that can be used by internals solvers. Given a list in `.nmr` format, let us say the `1A57.nmr` file, the following example show how simple is to use this function.

## Example
```julia-repl
preprocessing("1A57.nmr")
```
"""
function preprocessing(file::String,opt="read")
	if opt == "read"
		nmrfile = readdlm("$(file)")
		I = Int64.(nmrfile[:,1])
		J = Int64.(nmrfile[:,2])
		n = last(I)
		vpath = generate_virtual_path(I,J)
		lenI = length(I)
		vadd = Vector{Vector{Int64}}(undef,I[lenI])
		for i=1:I[lenI]
			vadd[i]=[]
		end
		V = Vector{NMRInfo}(undef,lenI)
		i = 1
		while i<=lenI
			k = 3 
			while i+k<=lenI && I[i]==I[i+k] 
				push!(vadd[I[i]],J[i+k])
				k = k+1
			end
			i = i+k
		end

		for i=1:lenI
			V[i]=NMRInfo(nmrfile[i,5],nmrfile[i,7],nmrfile[i,8])
		end
		for i=1:lenI
			if I[i]!=J[i]
				push!(I,J[i])
				push!(J,I[i])
				push!(V,NMRInfo(V[i].dist,V[i].atom2,V[i].atom1))
			end
		end
		nmrt = NMRType(vpath,vadd,sparse(I,J,V),n)
	else
		error("Unidentified option or file")	
	end
	return nmrt
end

"""
```
generate_virtual_path :: Function
```
This function is an auxiliary function used to define a useful vector called in our context as virtual path. This vector allows to handle with re-order approach. 
"""
function generate_virtual_path(NMRdatavertex1::Vector{Int64},NMRdatavertex2::Vector{Int64})
	#	D = [NMRdata.vertex1 NMRdata.vertex2]
	virtual_path = [1,2,3,4]
	k = 5
	li = 4
	while k <= last(NMRdatavertex1)
		ind = findall(x->x==k,NMRdatavertex1)
		la = 0 
		if NMRdatavertex2[ind[1:3]] ==  virtual_path[li-2:li]
			la+=1
		else 
			if NMRdatavertex2[ind[1:2]] == virtual_path[li-2:li-1] 
				la+=2
				push!(virtual_path,NMRdatavertex2[ind[3]])
			elseif NMRdatavertex2[ind[1:2]] == virtual_path[li-1:li]
				la+=2
				push!(virtual_path,NMRdatavertex2[ind[3]])
			elseif NMRdatavertex2[ind[2:3]] == virtual_path[li-2:li-1]
				la+=2
				push!(virtual_path,NMRdatavertex2[ind[1]])
			elseif NMRdatavertex2[ind[2:3]] == virtual_path[li-2:li]
				la+=2
				push!(virtual_path,NMRdatavertex2[ind[1]])
			else
				la+=4
				append!(virtual_path,NMRdatavertex2[ind[1:3]])
			end

		end
		li = li+la
		push!(virtual_path,k)
		k += 1
	end
	return virtual_path
end




# these types are mandatory indepedent of used solver##############################
"""
```
ConformationSetup :: Type
```
This type is used to define the main options of the function used to conformation. Essentially, this type has three fields: precision, solver, allsolutions. 

precision ::  is associated to current internal solver which depends on some level do a pruning test.

solver :: is the name of function which define the solver. 

allsolution :: is a bool variable used to define the behavior of the method. Is set as true, the internal solver try to find all possible solutions. 

## Example

```julia-repl
options = ConformationSetup(0.00001,classicBP,true)
```
As return an element options was created an it will used in conformation function. It is a mandatory definition before run conformation function.

"""
struct ConformationSetup
	precision :: Float64
	virtual_precision :: Float64
	solver :: Function 
	evalLDE :: Bool
	allsolutions :: Bool
	function ConformationSetup(a,c,f,d=true,b=1.0e-6)
		return new(a,b,c,d,f)
	end
	function ConformationSetup(;precision=1.0e-6, virtual_precision=precision, solver, evalLDE=true, allsolutions=false)
		return new(precision,virtual_precision,solver,evalLDE,allsolutions)
	end
end

"""
```
AtomType :: Type
```
It is the most primitive type and it is used to store the spatial positions of atoms. It is a mutable type. This type has three field: .x, .y and .z
"""
mutable struct AtomType
	x :: Float64
	y :: Float64
	z :: Float64
	element :: String
	function AtomType(a,b,c,d="None")
		return new(a,b,c,d)
	end
end


"""
```
MoleculeType :: Type
```
It is used to store a (just one) solution find by conformation process. It is a mutable struct and has two fields: .atoms (used to store a vector of AtomType) and .lde (used to store the lde of an molecule)
"""
mutable struct MoleculeType
	atoms :: Vector{AtomType}
	lde :: Float64
end

function writefile(A::MoleculeType,file="SpatialPosition";format = "xyz")
	if isa(A,MoleculeType)
		B = [length(A.atoms)," "]
		for i=1:length(A.atoms)
			push!(B,[A.atoms[i].element,A.atoms[i].x,A.atoms[i].y,A.atoms[i].z])
		end
		open("$(file).$(format)",io) do io
			writedlm(io,B)
		end
	else
		error("Not recognize format")
	end
end

"""
```
Counter :: Type
```
It is a mutable type for counting operations in solvers. This can handle with counters for nodes, virtual positions, ddf, branchs and prunes. 
"""
mutable struct Counter
	node :: Vector{Int64}
	virtual_path :: Vector{Int64}
	ddf :: Vector{Int64}
	branch :: Int64
	prune :: Int64
end

"""
```
ConformationOutput :: Type
```
It is a mutable type used to store the output provided by conformation function. With this type is possible to handle with some importants elements given by the following fields:.number,.molecules,.elapsedtime,.bytes and .gctime
"""
mutable struct ConformationOutput
	solver :: Function
	number :: Int64
	molecules :: Dict{Int64,MoleculeType}
end


###################################################################################
# these functions are used by ClassicBP  ##########################################
"""
```
bondangle :: Function
```
This is an auxiliary function used by ClassicBP solver in order to compute the bond angle by cosine rule. As output the cosine and sine are given.
"""
function bondangle(d23,d24,d34)
	c = (-d24*d24 + d34*d34 + d23*d23)/(2.0*d34*d23)
	if c<-1.0
		return -1.0, 0
	elseif c>1.0
		return 1.0, 0.0
	end
	return c,sqrt(1.0-c^2)
end

"""
```
badtorsionangle :: Function
```
This is an auxiliary function used by classicBP solver in order to compute the torsion angle. As output the cosine and sine of the torsion angle are given.
"""
function badtorsionangle(d12,d13,d14,d23,d24,d34)#i=4,...,n
	b = (d24*d24 + d23*d23 - d34*d34)/(2.0*d24*d23)
	c = (d12*d12 + d23*d23 - d13*d13)/(2.0*d12*d23)
	if (abs(b) > 1.0  || abs(c) > 1.0) 
		return -2
	end
	valc = ((d12*d12 + d24*d24 - d14*d14)/(2.0*d12*d24) - b*c)/(sqrt((1.0 - b*b)*(1.0 - c*c)))
	if (valc < -1.0)  
		return -1.0, 0.0
	elseif (valc >  1.0)  
		return 1.0, 0
	end
	return valc,sqrt(1.0-valc*valc)
end


"""
```
torsionmatrix :: Function
```
This is an auxiliary function  used by classicBp solver to compute the torsion array.
"""
function torsionmatrix(cosθ,sinθ,cosω,sinω,d34,B,sign::Bool)
	if sign == true

		B=zeros(4,4)
		B[1,1] = -cosθ
		B[1,2] = -sinθ
		B[1,4] = -d34*cosθ
		B[2,1] = sinθ*cosω
		B[2,2] = -cosθ*cosω
		B[2,3] = -sinω
		B[2,4] = d34*B[2,1]
		B[3,1] = sinθ*sinω
		B[3,2] = -cosθ*sinω
		B[3,3] = cosω
		B[3,4] = d34*B[3,1] 
		B[4,4] = 1
	else

		B[2,3] = -B[2,3]
		B[3,1] = -B[3,1]
		B[3,2] = -B[3,2]
		B[3,4] = -B[3,4]	
	end
	return B
end

function torsionmatrix(cosθ,sinθ,cosω,sinω,d34)
	B=zeros(4,4)
	B[1,1] = -cosθ
	B[1,2] = -sinθ
	B[1,4] = -d34*cosθ
	B[2,1] = sinθ*cosω
	B[2,2] = -cosθ*cosω
	B[2,3] = -sinω
	B[2,4] = d34*B[2,1]
	B[3,1] = sinθ*sinω
	B[3,2] = -cosθ*sinω
	B[3,3] = cosω
	B[3,4] = d34*B[3,1] 
	B[4,4] = 1
	return B
end
function torsionmatrix(C)
	B = copy(C)
	B[2,3] = -B[2,3]
	B[3,1] = -B[3,1]
	B[3,2] = -B[3,2]
	B[3,4] = -B[3,4]	
	return B
end
"""
```
pruningtest :: Function
```
This functions is an auxiliary function used to test if some molecule is 
feasible or not.
"""
function pruningtest(v::MoleculeType,i::Int64,D::NMRType,ε::Float64) :: Bool
	if isempty(D.additional_distance[i])

		@debug "result of prunning true"
		return true
	else
		for j in D.additional_distance[i]
			dij =  (v.atoms[i].x-v.atoms[j].x)^2+(v.atoms[i].y-v.atoms[j].y)^2 +(v.atoms[i].z-v.atoms[j].z)^2
			if (D.info[i,j].dist^2 -dij)^2 >ε
				@debug "result of prunning false"
				return false
			end
		end
	end

	@debug "result of prunning true"
	return true
end


####################################################################################
# these functions can be used by any solver#########################################
"""
```
LDE :: Function
```
This function is an auxiliary function used to compute the LDE of some molecule.
This function change the .lde field of the molecule.
"""
function LDE(v::MoleculeType,D::NMRType)
	dij=0.0
	nne = findnz(D.info)
	num_nne = length(nne[1])
	maxε = 0.0
	for k=1:num_nne
		i=nne[1][k]
		j=nne[2][k]
		aux = (D.info[i,j].dist^2-((v.atoms[i].x-v.atoms[j].x)^2+(v.atoms[i].y-v.atoms[j].y)^2 +(v.atoms[i].z-v.atoms[j].z)^2))^2
		maxε = max(maxε,aux)
		dij = dij+(aux)/D.info[i,j].dist
	end
	#	println("max ε = $(maxε)")
	#	if nad>0.0
	v.lde = dij/(2.0*num_nne)
	#	else
	#		return dij
	#	end
end

"""
```
outputfilter
```
This function is especially useful when we need to filter an information about a solution in ConformationOutput. For example, let us suppose that we want to compute the worst LDE measure. So, we need to type:

julia> a = conformation(data,option)

julia> outputfilter(a,"lde")
"""
function outputfilter(a::ConformationOutput, option = "lde")
	mol = a.molecules
	num = a.number
	if option == "lde"
		vlde = zeros(num)

		for i=1:num
			vlde[i]=mol[i].lde
		end
		return minimum(vlde)
	end
	if option == "xyz"
		nl = length(a.molecules[1].atoms)
		A=Array{Vector{Float64},2}(undef,nl,num)
		for j=1:num
			for i=1:nl
				A[i,j]=[a.molecules[j].atoms[i].x, a.molecules[j].atoms[i].y,
						a.molecules[j].atoms[i].z]
			end
		end
		return A
	end
end

"""
```
build_distance_matrix
```
This function can be useful to visualize a the distance array associated to an
MoleculeType solution

## Example

```julia-repl
v = conformation(data,opt)
build_distance_array(v.molecules[1].atoms)

```
"""
function build_distance_matrix(v::Array{AtomType,1})
	n=length(v)
	D=zeros(n,n)
	for i=1:n
		for j=1:n
			D[i,j] = sqrt((v[i].x-v[j].x)^2+(v[i].y-v[j].y)^2+(v[i].z-v[j].z)^2)
		end
	end
	return D
end

function Base.copy(A::AtomType)
	return AtomType(A.x,A.y,A.z,A.element)
end

function Base.copy(M::MoleculeType)
	nl = length(M.atoms)
	vm = Vector{AtomType}(undef,nl)
	for i=1:nl
		vm[i]=copy(M.atoms[i])
	end
	return MoleculeType(vm,M.lde)
end

function Base.:≈(A::MoleculeType,B::MoleculeType)
	na = length(A.atoms)
	nb = length(B.atoms)
	v = zeros(3)
	if na == nb 
		for i=1:na
			v = [A.atoms[i].x-B.atoms[i].x,A.atoms[i].y-B.atoms[i].y,A.atoms[i].z-B.atoms[i].z]
			if norm(v)>1.0e-4
				return false
			end
		end
	end
	return true
end
####################################################################################
function symsparse(I,J,v)
	for i=1:length(I)
		if I[i]!=J[i]
			push!(I,J[i])

			push!(J,I[i])
			push!(v,v[i])
		end
	end
	return sparse(I,J,v)
end

# these functions are used by quaternion_bp#########################################
function qbondangle(d23,d24,d34)
	c = (-d24*d24+ d34*d34+ d23*d23)/(4.0*d34*d23)
	if c<-0.5
		return 0.0, 1.0
	elseif c>0.5
		return 1.0, 0.0
	end
	return sqrt(0.5 + c),sqrt(0.5 - c)
end

function qtorsionangle(d12,d13,d14,d23,d24,d34)
	d12m = d12*d12 
	d23m = d23*d23
	d24m = d24*d24
	d34m = d34*d34
	dtil3 = d12m + d23m - d13*d13
	dtil2 = d24m + d23m - d34m        
	valc = (d23m*(d12m + d24m - d14*d14) - 0.5*dtil3*dtil2)/sqrt((4.0*d12m*d23m - dtil3*dtil3)*(4.0*d23m*d24m - dtil2*dtil2))
	if (valc < -0.5)  
		return 0.0, 1.0
	elseif (valc >  0.5)  
		return 1.0, 0.0
	end
	return sqrt(0.5+valc),sqrt(0.5-valc)
end

# Quaternion small library
mutable struct Quaternion
	s :: Float64
	v1 :: Float64
	v2 :: Float64
	v3 :: Float64
end

function Base.copy(Q::Quaternion)
	return Quaternion(Q.s, Q.v1, Q.v2, Q.v3)
end

#function qsign(q::Quaternion)
#	return Quaternion(-q.s, -q.v1, -q.v2, -q.v3)
#end

#function qsum(q::Quaternion,w::Quaternion)
#	return Quaternion(q.s + w.s, q.v1 + w.v1, q.v2 + w.v2, q.v3 + w.v3)
#end

#function qminus(q::Quaternion,w::Quaternion)
#	return Quaternion(q.s - w.s, q.v1 - w.v1, q.v2 - w.v2, q.v3 - w.v3)
#end

function qprod(q::Quaternion,w::Quaternion)
	return  Quaternion(q.s * w.s - q.v1 * w.v1 - q.v2 * w.v2 - q.v3 * w.v3,
					   q.s * w.v1 + q.v1 * w.s + q.v2 * w.v3 - q.v3 * w.v2,
					   q.s * w.v2 - q.v1 * w.v3 + q.v2 * w.s + q.v3 * w.v1,
					   q.s * w.v3 + q.v1 * w.v2 - q.v2 * w.v1 + q.v3 * w.s)
end

function qprod(q::Quaternion,a::Float64,b::Float64,c::Float64,d::Float64)
	return  Quaternion(q.s * a - q.v1 * b - q.v2 * c - q.v3 * d,
					   q.s * b + q.v1 * a + q.v2 * d - q.v3 * c,
					   q.s * c - q.v1 * d + q.v2 * a + q.v3 * b,
					   q.s * d + q.v1 * c - q.v2 * b + q.v3 * a)
end

#function conj(q::Quaternion)
#	return  Quaternion(q.s, -q.v1, -q.v2, -q.v3)
#end



function rotopt(Q::Quaternion,t::Float64)
	sl = 2.0*t
	return Quaternion(0.0, sl*(Q.s*Q.s + Q.v1*Q.v1 -0.5) ,sl*(Q.v2*Q.v1 + Q.v3*Q.s), sl*(Q.v3*Q.v1 - Q.v2*Q.s))
end

# to be fair with memory acess in comparations
function prodmatrix(A::Array{Float64,2},B::Array{Float64,2})
	C=zeros(4,4)
	C[4,4] = 1.0
	for i=1:3
		C[i,1] = A[i,1]*B[1,1] + A[i,2]*B[2,1] + A[i,3]*B[3,1] 
		C[i,2] = A[i,1]*B[1,2] + A[i,2]*B[2,2] + A[i,3]*B[3,2] 
		C[i,3] = A[i,2]*B[2,3] + A[i,3]*B[3,3]
		C[i,4] = A[i,1]*B[1,4] + A[i,2]*B[2,4] + A[i,3]*B[3,4] + A[i,4]
	end
	return C 
end
#[24,33,0,0]
#[48,64,0,0]
"""
	convert_to_dataframe

A function used to convert an array of AtomType elements in a dataframe

# Examples
```
julia-repl
julia> sol = conformation(data,options) # assuming all solutions founded

julia> convert_to_dataframe(sol.molecules[1].atoms)

returns a dataframe object. 
```
"""
function convert_to_dataframe(A::Array{AtomType,1})
	n = length(A)
	df = DataFrame()
	df.atom = []
	df.x = []
	df.y = []
	df.z = []
	for k = 1:n
		push!(df.atom,A[k].element)
		push!(df.x,A[k].x)
		push!(df.y,A[k].y)
		push!(df.z,A[k].z)
	end
	return df
end

