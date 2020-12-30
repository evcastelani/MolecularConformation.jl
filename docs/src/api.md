# API

We can summarize the main functions of this package as follows.

## Types

```@docs
NMRInfo
NMRType
AtomType
MoleculeType
ConformationSetup
ConformationOutput
```


## Main functions

```@docs
preprocessing(file::String,opt="read")
conformation(NMRdata::NMRType,cs::ConformationSetup)
classicBP(NMRdata :: NMRType,ε :: Float64,virtual_ε :: Float64,allmol :: Bool)
LDE(v::MoleculeType,D::NMRType)
build_distance_matrix(v::Array{AtomType,1})
convert_to_dataframe(A::Array{AtomType,1})
```



## Secundary functions


This functions are used as auxiliary functions by internal solvers.

```@docs
generate_virtual_path(NMRdatavertex1::Vector{Int64},NMRdatavertex2::Vector{Int64})
bondangle(d23,d24,d34)
torsionangle(d12,d13,d14,d23,d24,d34)
badtorsionangle(d12,d13,d14,d23,d24,d34)
torsionmatrix(cosθ,sinθ,cosω,sinω,d34,B,sign::Bool)
pruningtest(v::MoleculeType,i::Int64,D::NMRType,ε::Float64,count::Vector{Int64})
outputfilter(a::ConformationOutput, option = "lde")

```

