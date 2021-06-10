using MolecularConformation, DelimitedFiles,BenchmarkTools, Dates

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10

include("rmsd.jl")
"""
The perform function is used to run tests in order to compare algorithms. 
There are some dependences to use this functions like BenchmarkTools and
Dates. 

## Example
```julia-repl
julia> include("perform2.jl")

julia> perform(Second(3.0))
```
Note that now a timeout can be used. Some informations are displayed in 
order to show how things are going on. 
When the execution is finished, three files are build: not_solved.csv,
table_general.tex and table_operations.tex. By default these files are 
write because the option opwrite is defined as "a", that is, append. 
If you want redefine the files, a new definition need to be set. 

## More examples
```julia-repl
julia> include("perform2.jl")

julia> perform(Second(3.0),"w")

julia> perform(Minute(2)) # to timeout of 2 minutes
```

"""
function perform(limit_time,opwrite::String="a",f::Function=median)
	# setup 
	optq = ConformationSetup(0.000001,quaternionBP,false)
	optc = ConformationSetup(0.000001,classicBP,false)
	data = preprocessing("toyinstance.nmr")
	# first run (to optimize)
	solq = conformation(data,optq,limit_time)
	solc = conformation(data,optc,limit_time)
	list_of_problems = ["pdb1a57", "pdb1b4c", "pdb1ba5"]

	io = open("not_solved.csv",opwrite)
	iogen = open("table_general.tex",opwrite)
	ioop  = open("table_operations.tex",opwrite)

	println(iogen,"Problem & Method & LDE & PT & Num. Sol & rmsd & Improv \\\\")
	println(ioop, "Problem & Method & op. nodes & op. virtual path & op. ddf & Num. Branch & Num. Prune & Improv. nop & Improv. total \\\\")

	c = 1e-9
	for prob in list_of_problems
		# before make the benchmark we need to verify if is possible solve the problem
		# considering the limit_time
		data = preprocessing(string(prob,".nmr"))
		bench = true
		try 
			solq = conformation(data,optq,limit_time)
		catch
			bench = false
		end
		if bench == true
			print(" 🔔 The problem $(prob) can be solved within the timeout. \n") 
			print(" 👏 The benchmark will be performed for both methods! \n")
			print(" 🕐 Running benchmark using QuaternionBP... \n")
			bchq = @benchmark conformation($(data),$(optq),$(limit_time)+Second(10.0))
			print(" 🏁 The benchmark using QuaternionBP in $(prob) was done!\n")
			LDEq = outputfilter(solq,"lde")
			PTq = f(bchq).time*c

			print(" 🕞 Running benchmark using ClassicBP... \n")
			bchc = @benchmark conformation($(data),$(optc),$(limit_time)+Second(10.0))
			print(" 🏁 The benchmark using ClassicBP in $(prob) was done!\n")
			print(" \n")

			solq = conformation(data,optq,limit_time+Second(10.0))
			solc = conformation(data,optc,limit_time+Second(10.0))	
			LDEc = outputfilter(solc,"lde")
			PTc = f(bchc).time*c
			improv = (-1.0+PTc/PTq)*100
			rmsd = evalrmsd(solc,solq)[2]
		
			println(iogen,"$(prob) & quaternionBP & $(LDEq) & $(PTq) & $(solq.number) & $(rmsd) & $(improv)\\\\")
			println(iogen,"$(prob) & classicBP & $(LDEc) & $(PTc) & $(solc.number) & - & - \\\\")

			improv_nop = (-1.0 .+ solc.nop.node./solq.nop.node).*100
			improv_total = (-1.0 + sum(solc.nop.node)/sum(solq.nop.node))*100
			println(ioop,"$(prob) & quaternionBP & $(solq.nop.node) & $(solq.nop.virtual_path) & $(solq.nop.ddf) & $(solq.nop.branch) & $(solq.nop.prune) & $(improv_nop) & $(improv_total) \\\\")
			println(ioop,"$(prob) & classicBP & $(solc.nop.node) & $(solc.nop.virtual_path) & $(solc.nop.ddf) & $(solc.nop.branch) & $(solc.nop.prune) & - & - \\\\")
			
		else
			print(" 💥 The problem $(prob) cannot be solved within the timeout.\n")
			print(" \n")
			println(io, prob )
		end

	end
	close(io)
	close(iogen)
	close(ioop)


end

function vec2array(vec)
	len = length(vec)
	A=zeros(len,3)
	for i=1:len
		A[i,1] = vec[i][1]
		A[i,2] = vec[i][2]
		A[i,3] = vec[i][3]
	end
	return A
end


function evalrmsd(sols1::ConformationOutput,sols2::ConformationOutput)
	coordsols1 = outputfilter(sols1,"xyz")
	coordsols2 = outputfilter(sols2,"xyz")
	valrmsd = [1000.0,1000.0]
	#display(coordfile)
	for k=1:sols1.number
		coords1 = vec2array(coordsols1[:,k])
		coords2 = vec2array(coordsols2[:,k])
		valrmsdaux = rmsd(coords1,coords2)
		if valrmsdaux[2] <= valrmsd[2]
			valrmsd = [valrmsdaux[1],valrmsdaux[2]]
		end
	end
	return valrmsd
end
