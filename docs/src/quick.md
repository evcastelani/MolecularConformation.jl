## Installation

Currently, `MolecularConformation.jl` is not an registered package in Julia database, consequently, the best way to install is type `]` in Julia REPL and then

```
pkg> add "https://github.com/evcastelani/MolecularConformation.jl"
```

If you want to setup some especific branch, is necessary just include the term `#nameofbranch` in the end of previous string. 

## Basic Usage


First of all, we need to load the package. 

```@repl tutorial
	using MolecularConformation
```

In order to conformalize some structure,  we need to get a list in `.nmr` format. This kind of list can get, for example, [here](https://github.com/mucherino/mdjeep/tree/master/instances/0.3/proteinSet1). 

!!! warning " 'Warning' "
    Currently, `MolecularConformation.jl` works only with file in `.nmr ` format. 


Since we have the list in this format, we must transform the data to `NMRType` format. The simplest way to do this is through the `preprocessing` function. An `NMRType` element contains information about order, additional distances and other elements that can be used to optimize the solver. The `preprocessing` function can handle with the `.nmr` file and provides a `NMRType` with all necessary information.

```@repl tutorial
	data_example = preprocessing("toyinstance.nmr")
```

Note that the type of `data_example` is a `NMRType`

```@repl tutorial
	typeof(data_example)
```

which means that `data_example` is has four fields.

```@repl tutorial
	fieldnames(NMRType)
```

One of most important fields is `.info`. This field is a symmetric `CSVSparseMatrix` where each input value is a `NMRInfo` element. This kind of structure provides a fast acess and considerable storage savings. A `NMRInfo` element contains information about the `distance` and two bond atoms.  

```@repl tutorial
	fieldnames(NMRInfo)
```
 
Once the type of data this package handles is well established, let's look at the rest of the settings until the main function is executed. The next step is to define the options to be used for the solver. For this, we need to create a `ConformationSetup`. This kind of element has five fields but two are optional. The field `virtual_precision` is one of optional argument and it is used to validated a step of use of re-order strategy. By default `virtual_precision` is set like `precision` field. 

```@repl tutorial
	fieldnames(ConformationSetup)
```
The `precision`  will be used by the pruning test of the solver. The `solver` is the name of the function that was choose and `allsolution` is a bolean parameter for detect all solution if it is set as true. The `evalLDE` is another optional argument and 

```@repl tutorial
	opt = ConformationSetup(0.000001,MolecularConformation.classicBP,false)
```
  
Now all that is required to execute the algorithm is ready and therefore, let us run the `conformation` function. 

```@repl tutorial
	sol = conformation(data_example,opt)
```

In our example, the solution is storage in `sol` whose type is `ConformationOutput`.

```@repl tutorial
	typeof(sol)
```

The `ConformationOutput` type, has the following fields.

```@repl tutorial
	fieldnames(ConformationOutput)
```

The field `solver` is the used solver defined in setup. The field `number` corresponds to the number of solutions. The field `molecules` is a dictionary of where associated to each key we have a solution of `MoleculeType`.

```@repl tutorial
	typeof(sol.molecules[1])
```

A `MoleculeType` has two fields.

```@repl tutorial
	fieldnames(MoleculeType)
```

In our example this has the following values

```@repl tutorial
	sol.molecules[1].atoms
	sol.molecules[1].lde
```

The last field of `ConformationOutput` is called `nop` which means `number of operations` and it is self-explanatory and used to just to analyze the performance of the method.

## About the solvers

Currently, just three solvers are available: `classicBP`, `classicBPOpt` and `quaternionBP`.
