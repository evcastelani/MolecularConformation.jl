# these types are mandatory indepedent of used solver##############################
#
"""
```
NMRinfo
```
This type contains the information for NMRtype. In this point of the project it contains two fields: dist and typeatom. 
"""
struct NMRInfo
	dist::Float64
	typeatom::String
end

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
	virtual_path :: Vector{Int64}
	additional_distance::Vector{Vector{Int64}}
	info :: SparseMatrixCSC{NMRInfo,Int64}
	dim :: Int64
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
			V[i]=NMRInfo(nmrfile[i,5],nmrfile[i,7])
		end
		for i=1:lenI
               		if I[i]!=J[i]
                   		push!(I,J[i])
                        	push!(J,I[i])
                  	 	push!(V,V[i])
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
generate_virtual_path
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


###################################################################################
# these functions are used by bp_classical##########################################
function bondangle(d23,d24,d34)#i=3,...,n
	c = (-d24^2 + d34^2 + d23^2)/(2.0*d34*d23)
#	c = (-D[i-2,i]^2+ D[i-1,i]^2+ D[i-2,i-1]^2)/(2.0*D[i-1,i]*D[i-2,i-1])
	if c<-1.0
		c=-1.0
	end
	if c>1.0
		c=1.0
	end
	s = sqrt(1.0-c^2)
	return c,s
end

function torsionangle(d12,d13,d14,d23,d24,d34)#i=4,...,n
	#d12=D[i-3,i-2]
	#d13=D[i-3,i-1]
	#d14=D[i-3,i]
	#d23=D[i-2,i-1]
	#d24=D[i-2,i]
	#d34=D[i-1,i]
	a = d12*d12 + d24*d24 - d14*d14
	a = a/(2.0*d12*d24)
	b = d24*d24 + d23*d23 - d34*d34        
	b = b/(2.0*d24*d23)
	c = d12*d12 + d23*d23 - d13*d13
	c = c/(2.0*d12*d23)
	e = 1.0 - b^2;
	f = 1.0 - c^2;
	if (e < 0.0 || f < 0.0)  
		@debug "some problem in torsion angle"
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

function torsionmatrix(cosθ,sinθ,cosω,sinω,d34,sign::Char)
	if sign == '+'
	
		B=zeros(4,4)
		B[1,1] = -cosθ
		B[1,2] = -sinθ
		B[1,4] = -d34*cosθ
		B[2,1] = sinθ*cosω
		B[2,2] = -cosθ*cosω
		B[2,3] = -sinω
		B[2,4] = d34*sinθ*cosω
		B[3,1] = sinθ*sinω
		B[3,2] = -cosθ*sinω
		B[3,3] = cosω
		B[3,4] = d34*sinθ*sinω 
		B[4,4] = 1
	else
	
		B=zeros(4,4)
		B[1,1] = -cosθ
		B[1,2] = -sinθ
		B[1,4] = -d34*cosθ
		B[2,1] = sinθ*cosω
		B[2,2] = -cosθ*cosω
		B[2,3] = sinω
		B[2,4] = d34*sinθ*cosω
		B[3,1] = -sinθ*sinω
		B[3,2] = cosθ*sinω
		B[3,3] = cosω
		B[3,4] = -d34*sinθ*sinω 
		B[4,4] = 1
	
	end
	return B
end

"""
```
pruningtest
```
This functions is an auxiliary function used to test if some molecule is 
feasible or not.
"""
function pruningtest(v::MoleculeType,i::Int,D::NMRType,ε::Float64)
	if isempty(D.additional_distance[i])
		return 1
	else
		for j in D.additional_distance[i]
			dij =  (v.atoms[i].x-v.atoms[j].x)^2+(v.atoms[i].y-v.atoms[j].y)^2 +(v.atoms[i].z-v.atoms[j].z)^2
			if (D.info[i,j].dist -dij)^2 >ε
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
This function change the .lde field of the molecule.
"""
function LDE(v::MoleculeType,D::NMRType)
	dij=0.0
	nne = findnz(D.info)
	num_nne = length(nne[1])
	for k=1:num_nne
		i=nne[1][k]
		j=nne[2][k]
		dij = dij+(D.info[i,j].dist^2-((v.atoms[i].x-v.atoms[j].x)^2+(v.atoms[i].y-v.atoms[j].y)^2 +(v.atoms[i].z-v.atoms[j].z)^2))^2/D.info[i,j].dist
	end
	
#	if nad>0.0
	return dij/(2.0*num_nne)
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


