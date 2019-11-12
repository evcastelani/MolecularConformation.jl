# these types are mandatory indepedent of used solver##############################
"""
```
NMRType
```
This type is used as input data. The present type was motivated by https://github.com/mucherino/mdjeep and it is used to avoid the definition of distance matrix. Once the .nmr file is read, this type is produced as output. Therefore, the simplest construction of this type is done using the nmr function as in the following example.

## Example

```julia-repl
nmr("1a57")
```
As return an element of NMRtype is created. We are assuming that a .nmr file as in  https://github.com/mucherino/mdjeep is given.
"""
struct NMRType
	vertex1 :: Vector{Int64}
	vertex2 :: Vector{Int64}
	gvertex1 :: Vector{Int64}
	gvertex2 :: Vector{Int64}
	lowerbound :: Vector{Float64}
	upperbound :: Vector{Float64}
	label1 :: Vector{String}
	label2 :: Vector{String}
	labelg1 :: Vector{String}
	labelg2 :: Vector{String}
end

"""
```
nmr
```
It is a function used to read a PBD file in format .nmr or .mdjeep. Just one option is avaiable: read.
"""
function nmr(file::String,opt="read")
	if opt == "read"
		nmrfile = readdlm("$(file)")
		nmrt = NMRType(nmrfile[:,1],nmrfile[:,2],nmrfile[:,3],nmrfile[:,4],nmrfile[:,5],nmrfile[:,6],nmrfile[:,7],nmrfile[:,8],nmrfile[:,9],nmrfile[:,10])
	else
		error("Unidentified option or file")	
	end
	return nmrt
end


# these types are mandatory indepedent of used solver##############################
"""
```
ConformationSetup
```
This type is used to define the main options of the function used to conformation.Essentially, this type has four fields: precision,cutoff,solver,allsolutions, reorder, interval

## Example

```julia-repl
options = ConformationSetup(0.001,5.5,classical_bp,true)
```
As return an element options was created an it will used in conformation function. It is a mandatory definition before run conformation function.

"""
mutable struct ConformationSetup
	precision :: Float64
	cutoff :: Float64
	solver :: Function 
	allsolutions :: Bool
	reorder :: Bool
	interval :: Bool
end

"""
```
AtomType
```
It is the most primitive type and it is used to store the spatial positions of atoms. It is a mutable type. This type has three field: .x, .y and .z
"""
mutable struct AtomType
	x :: Float64
	y :: Float64
	z :: Float64
end


"""
```
MoleculeType
```
It is used to store a (just one) solution find by conformation process. It is a mutable struct and has two fields: .atoms (used to store a vector of AtomType) and .lde (used to store the lde of an molecule)
"""
mutable struct MoleculeType
	atoms :: Vector{AtomType}
	lde :: Float64
end

"""
```
ConformationOutput
```
It is a mutable type used to store the output provided by conformation function. With this type is possible to handle with some importants elements given by the following fields:.number,.molecules,.elapsedtime,.bytes and .gctime
"""
mutable struct ConformationOutput
	number :: Int64
	molecules :: Dict{Int64,MoleculeType}
	elapsedtime :: Any
	bytes :: Any 
	gctime :: Any
end

"""
```
generate_virtual_path
```
This function is an auxiliary function used to define a useful vector called in our context as virtual path. This vector allows to handle with re-order approach. 
"""
function generate_virtual_path(NMRdata::NMRType)
#	D = [NMRdata.vertex1 NMRdata.vertex2]
	virtual_path = [1,2,3,4]
	k = 5
	li = 4
	while k <= last(NMRdata.vertex1)
		ind = findall(x->x==k,NMRdata.vertex1)
		la = 0 
		if NMRdata.vertex2[ind[1:3]] ==  virtual_path[li-2:li]
			la+=1
		else 
			if NMRdata.vertex2[ind[1:2]] == virtual_path[li-2:li-1] 
				la+=2
				push!(virtual_path,NMRdata.vertex2[ind[3]])
			elseif NMRdata.vertex2[ind[1:2]] == virtual_path[li-1:li]
				la+=2
				push!(virtual_path,NMRdata.vertex2[ind[3]])
			elseif NMRdata.vertex2[ind[2:3]] == virtual_path[li-2:li-1]
				la+=2
				push!(virtual_path,NMRdata.vertex2[ind[1]])
			elseif NMRdata.vertex2[ind[2:3]] == virtual_path[li-2:li]
				la+=2
				push!(virtual_path,NMRdata.vertex2[ind[1]])
			else
				la+=4
				append!(virtual_path,NMRdata.vertex2[ind[1:3]])
			end

		end
		li = li+la
		push!(virtual_path,k)
		k += 1
	end
	return virtual_path
end



###################################################################################
# these functions are used by bp_classical##########################################
function bondangle(i,D)#i=3,...,n
	c = (-D[i-2,i]^2+ D[i-1,i]^2+ D[i-2,i-1]^2)/(2.0*D[i-1,i]*D[i-2,i-1])
	if c<-1.0
		c=-1.0
	end
	if c>1.0
		c=1.0
	end
	s = sqrt(1.0-c^2)
	return c,s
end

function torsionangle(i,D)#i=4,...,n
	d12=D[i-3,i-2]
	d13=D[i-3,i-1]
	d14=D[i-3,i]
	d23=D[i-2,i-1]
	d24=D[i-2,i]
	d34=D[i-1,i]
	a = d12*d12 + d24*d24 - d14*d14
	a = a/(2.0*d12*d24)
	b = d24*d24 + d23*d23 - d34*d34        
	b = b / (2.0*d24*d23)
	c = d12*d12 + d23*d23 - d13*d13
	c = c / (2.0*d12*d23)
	e = 1.0 - b^2;
	f = 1.0 - c^2;
	if (e < 0.0 || f < 0.0)  
		return -2
	end
	e = sqrt(e)
	f = sqrt(f)
	valc = (a - b*c)/(e*f)
	if (valc < -1.0)  
		valc = -1.0
	end
	if (valc >  1.0)  
		valc =  1.0
	end
	vals=sqrt(1.0-valc^2)
	return valc,vals
end

function torsionmatrix(i,D,sign::Char)
	if sign == '+'
		ba = bondangle(i,D)
		ta = torsionangle(i,D)
	
		B=zeros(4,4)
		B[1,1] = -ba.cosθ
		B[1,2] = -ba.sinθ
		B[1,4] = -D[i,i-1]*ba.cosθ
		B[2,1] = ba.sinθ*ta.cosω
		B[2,2] = -ba.cosθ*ta.cosω
		B[2,3] = -ta.sinω
		B[2,4] = D[i,i-1]*ba.sinθ*ta.cosω
		B[3,1] = ba.sinθ*ta.sinω
		B[3,2] = -ba.cosθ*ta.sinω
		B[3,3] = ta.cosω
		B[3,4] = D[i,i-1]*ba.sinθ*ta.sinω 
		B[4,4] = 1
	else
		ba = bondangle(i,D)
		ta = torsionangle(i,D)
	
		B=zeros(4,4)
		B[1,1] = -ba.cosθ
		B[1,2] = -ba.sinθ
		B[1,4] = -D[i,i-1]*ba.cosθ
		B[2,1] = ba.sinθ*ta.cosω
		B[2,2] = -ba.cosθ*ta.cosω
		B[2,3] = ta.sinω
		B[2,4] = D[i,i-1]*ba.sinθ*ta.cosω
		B[3,1] = -ba.sinθ*ta.sinω
		B[3,2] = ba.cosθ*ta.sinω
		B[3,3] = ta.cosω
		B[3,4] = -D[i,i-1]*ba.sinθ*ta.sinω 
		B[4,4] = 1
	
	end
	return B
end

"""
```
pruningteste
```
This functions is an auxiliary function used to test if some molecule is 
feasible or not.
"""
function pruningtest(v::MoleculeType,i::Int,D::Array{Float64,2},ε::Float64,ndiag::Int)
	if ndiag == 0
		initj = 1
	else
		initj = max(1,i-ndiag)
	end
	for j=initj:i-1
		if D[i,j]>0.0
			dij = (v.atoms[i].x-v.atoms[j].x)^2+(v.atoms[i].y-v.atoms[j].y)^2 +(v.atoms[i].z-v.atoms[j].z)^2
			if (D[i,j]^2-dij)^2 > ε
				return 0 
			end
		end
	end
	return 1 
end

#unction perfpruningtest(v::MoleculeType,i::Int,D::Array{Float64,2},ε::Float64,ndiag::Int)
#       
#       for j=max(1,i-ndiag):i-1
#       	if D[i,j]>0.0
#       		dij = (v.atoms[i].x-v.atoms[j].x)^2+(v.atoms[i].y-v.atoms[j].y)^2 +(v.atoms[i].z-v.atoms[j].z)^2
#       		if (D[i,j]^2-dij)^2 > ε
#       			return 0 
#       		end
#       	end
#       end
#       return 1 
#nd

####################################################################################

# these functions can be used by any solver#########################################
"""
```
LDE
```
This function is an auxiliary function used to compute the LDE of some molecule.
This function change the .lde field of molecule.
"""
function LDE(v::MoleculeType,D::Array{Float64,2},n::Int,nad::Int)
	dij=0.0
	for i=1:n
		for j=i+1:n
			if 0.0<D[i,j]
				dij = dij+(D[i,j]^2-((v.atoms[i].x-v.atoms[j].x)^2+(v.atoms[i].y-v.atoms[j].y)^2 +(v.atoms[i].z-v.atoms[j].z)^2))^2/D[i,j]
			end
		end
	end
#	if nad>0.0
		return dij/nad
#	else
#		return dij
#	end
end

function Base.copy(A::AtomType)
	return AtomType(A.x,A.y,A.z)
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

# these functions are used by quaternion_bp#########################################
function qbondangle(i,D)#i=3,...,n
	c = (-D[i-2,i]^2+ D[i-1,i]^2+ D[i-2,i-1]^2)/(2.0*D[i-1,i]*D[i-2,i-1])
	if c<-1.0
		c=-1.0
	end
	if c>1.0
		c=1.0
	end
	cm = c/2.0 
	return sqrt(0.5 + cm),sqrt(0.5 - cm)
end

function qtorsionangle(i,D)#i=4,...,n
	d12=D[i-3,i-2]
	d13=D[i-3,i-1]
	d14=D[i-3,i]
	d23=D[i-2,i-1]
	d24=D[i-2,i]
	d34=D[i-1,i]
	a = d12*d12 + d24*d24 - d14*d14
	a = a/(2.0*d12*d24)
	b = d24*d24 + d23*d23 - d34*d34        
	b = b / (2.0*d24*d23)
	c = d12*d12 + d23*d23 - d13*d13
	c = c / (2.0*d12*d23)
	e = 1.0 - b^2;
	f = 1.0 - c^2;
	if (e < 0.0 || f < 0.0)  
		return -2
	end
	e = sqrt(e)
	f = sqrt(f)
	valc = (a - b*c)/(e*f)
	if (valc < -1.0)  
		valc = -1.0
	end
	if (valc >  1.0)  
		valc =  1.0
	end
	cm = valc/2.0
	return sqrt(0.5+cm),sqrt(0.5-cm)
end


function prodmatrix(A::Array{Float64,2},B::Array{Float64,2})
	C=zeros(4,4)
	for i=1:4
		for j=1:4
			for k=1:4
				C[i,j]= C[i,j]+A[i,k]*B[k,j]
			end
		end
	end
	return C 
end

# Quaternion small library
mutable struct Quaternion
	s :: Float64
	v1 :: Float64
	v2 :: Float64
	v3 :: Float64
end

function qsign(q::Quaternion)
	return Quaternion(-q.s, -q.v1, -q.v2, -q.v3)
end

function qsum(q::Quaternion,w::Quaternion)
	return Quaternion(q.s + w.s, q.v1 + w.v1, q.v2 + w.v2, q.v3 + w.v3)
end

function qminus(q::Quaternion,w::Quaternion)
	return Quaternion(q.s - w.s, q.v1 - w.v1, q.v2 - w.v2, q.v3 - w.v3)
end
    
function qprod(q::Quaternion,w::Quaternion)
	return  Quaternion(q.s * w.s - q.v1 * w.v1 - q.v2 * w.v2 - q.v3 * w.v3,
                           q.s * w.v1 + q.v1 * w.s + q.v2 * w.v3 - q.v3 * w.v2,
                           q.s * w.v2 - q.v1 * w.v3 + q.v2 * w.s + q.v3 * w.v1,
                           q.s * w.v3 + q.v1 * w.v2 - q.v2 * w.v1 + q.v3 * w.s)
end

function conj(q::Quaternion)
	return  Quaternion(q.s, -q.v1, -q.v2, -q.v3)
end


function rot(Q::Quaternion,t::Float64)
	return qprod(Quaternion(-Q.v1*t,Q.s*t,Q.v3*t,-Q.v2*t),conj(Q))
end


