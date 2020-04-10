# Basic information to run `lavor instances`


`lavor instances` is a collection of artificial problems in order to test the robustness and efficiency of solvers implemented in `MolecularConformation.jl`. Consequently, this folder contains just the script to perform the algorithms and two class of problems given in folders `10` and `15`. There are more avaiable problems.  The reason to put the set in another place is to save storage in GitHub. 

Once the reason of this separation of folders is clear, the user can download more problems here. Put the folder related to some problem in the same place of the script `perform.jl` and `toyproblem.nmr`. To run just type:

```julia
julia> include("perform.jl")

julia> perform(10) # to problem given in folder 10

julia> perform(100) # to problem given in folder 100 and go on 
```
Our solvers are recursive. So, in some cases, can be necessary to extend the limit of memory in your system. For example, type in your terminal `limit -s 16384`. In general this number is enough. 

