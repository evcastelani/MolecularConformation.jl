using MolecularConformation, DelimitedFiles,BenchmarkTools, Dates, Printf, DelimitedFiles

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000

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
function perform(limit_time,opwrite::String="a",f::Function=median;list_of_problems="")
	# setup 
	optq = ConformationSetup(0.000001,quaternionBP,false)
	optc = ConformationSetup(0.000001,classicBP,false)
	data = preprocessing("toyinstance.nmr")
	# first run (to optimize)
	solq = conformation(data,optq,limit_time)
	solc = conformation(data,optc,limit_time)
	if length(list_of_problems) == 0 
		#list_of_problems = ["pdb1a57", "pdb1b4c", "pdb1ba5"]
		list_of_problems = ["pdb1a57", "pdb1b4c", "pdb1ba5", "pdb1d1n", "pdb1dp3", "pdb1du1", "pdb1eii", "pdb1fcl", "pdb1fd6", "pdb1hf9", "pdb1i2u", "pdb1i2v", "pdb1ijc", "pdb1jlz", "pdb1jw3", "pdb1k0v", "pdb1k2h", "pdb1k36", "pdb1k37", "pdb1kuw", "pdb1kz0", "pdb1kz2", "pdb1kz5", "pdb1lvz", "pdb1m4e", "pdb1ma2", "pdb1ma4", "pdb1ma5", "pdb1ma6", "pdb1mpe", "pdb1nd9", "pdb1ne5", "pdb1nmj", "pdb1o53", "pdb1oqp", "pdb1plw", "pdb1plx", "pdb1pv0", "pdb1pzr", "pdb1qlk", "pdb1r57", "pdb1ry3", "pdb1s4h", "pdb1s4j", "pdb1s6j", "pdb1sa8", "pdb1t2y", "pdb1t5q", "pdb1tot", "pdb1v6r", "pdb1v92", "pdb1vd7", "pdb1vd9", "pdb1vdb", "pdb1vpc", "pdb1wnk", "pdb1wo4", "pdb1wo5", "pdb1x60", "pdb1x9v", "pdb1y5c", "pdb1yx7", "pdb1yx8", "pdb1yxr", "pdb2a2y", "pdb2a4j", "pdb2adl", "pdb2adn", "pdb2ajj", "pdb2ajm", "pdb2ajn", "pdb2ajo", "pdb2akk", "pdb2bzb", "pdb2c0s", "pdb2dci", "pdb2eem", "pdb2fva", "pdb2fvf", "pdb2fxz", "pdb2g9j", "pdb2g9l", "pdb2h5m", "pdb2hep", "pdb2j0z", "pdb2j10", "pdb2j11", "pdb2jmy", "pdb2jn5", "pdb2jnk", "pdb2jpn", "pdb2jua", "pdb2jvd", "pdb2jwe", "pdb2jws", "pdb2jwu", "pdb2jxf", "pdb2jz5", "pdb2k2a", "pdb2k2f", "pdb2k36", "pdb2k37", "pdb2k3i", "pdb2k6s", "pdb2k7o", "pdb2kbm", "pdb2kdh", "pdb2kdl", "pdb2kdm", "pdb2kdp", "pdb2kdr", "pdb2kes", "pdb2kib", "pdb2kjn", "pdb2kjo", "pdb2kjr", "pdb2kko", "pdb2kl5", "pdb2klz", "pdb2ko1", "pdb2koz", "pdb2kp0", "pdb2ksg", "pdb2kt8", "pdb2kuh", "pdb2kwh", "pdb2kxa", "pdb2kyb", "pdb2l3m", "pdb2l3n", "pdb2l45", "pdb2l5r", "pdb2l6q", "pdb2l6r", "pdb2l98", "pdb2lci", "pdb2lde", "pdb2le2", "pdb2le7", "pdb2ler", "pdb2lgi", "pdb2lhc", "pdb2lhd", "pdb2lhe", "pdb2lhg", "pdb2lix", "pdb2lld", "pdb2lm9", "pdb2lmf", "pdb2ln3", "pdb2lo2", "pdb2lqp", "pdb2lr0", "pdb2lrh", "pdb2ls9", "pdb2lsa", "pdb2lss", "pdb2lt3", "pdb2lu6", "pdb2lwa", "pdb2lx0", "pdb2lxb", "pdb2lxy", "pdb2lxz", "pdb2lzf", "pdb2m1a", "pdb2m1j", "pdb2m2y", "pdb2m3f", "pdb2m4j", "pdb2m8m", "pdb2m8o", "pdb2m97", "pdb2m9g", "pdb2m9r", "pdb2mdv", "pdb2me1", "pdb2me2", "pdb2me3", "pdb2me4", "pdb2mg1", "pdb2mg2", "pdb2mg3", "pdb2mh8", "pdb2mhw", "pdb2mi1", "pdb2mi7", "pdb2mij", "pdb2mix", "pdb2mj1", "pdb2mj2", "pdb2mjg", "pdb2mji", "pdb2mle", "pdb2mlf", "pdb2mo5", "pdb2mpu", "pdb2mpz", "pdb2mra", "pdb2msu", "pdb2mty", "pdb2muh", "pdb2mvj", "pdb2mvt", "pdb2mvx", "pdb2mz6", "pdb2n00", "pdb2n0b", "pdb2n0c", "pdb2n0d", "pdb2n0e", "pdb2n0f", "pdb2n0g", "pdb2n0h", "pdb2n35", "pdb2n3a", "pdb2n41", "pdb2n4e", "pdb2n5q", "pdb2n67", "pdb2n6n", "pdb2n7i", "pdb2n7j", "pdb2n9b", "pdb2na9", "pdb2ndc", "pdb2nde", "pdb2ndk", "pdb2nvj", "pdb2p6j", "pdb2p81", "pdb2pqe", "pdb2rlg", "pdb2rlh", "pdb2rmy", "pdb2rnd", "pdb2roo", "pdb2rql", "pdb2rut", "pdb2ruv", "pdb2rux", "pdb2ruy", "pdb2rv3", "pdb2rv5", "pdb2y4q"]
	end

	io = open("not_solved.csv",opwrite)
    ioperf = open("results.csv",opwrite)
	iogen = open("table_general.tex",opwrite)
	
	println(iogen, "\\begin{xltabular}{\\textwidth}{r|rS[table-format=1.3e+2]S[table-format=1.4e+2]S[table-format=1.4e+2]S[table-format=-1.5]}
	\t\\caption{Results} \\label{tab:genResults}\\\\
	\t\\toprule
	\t\\multicolumn{1}{c}{\$\\begin{array}{c}
	\\text{Problem}\\\\ \\text{Size}
	\\end{array}\$} & Method & \\multicolumn{1}{c}{LDE} & \\multicolumn{1}{c}{RMSD} & \\multicolumn{1}{c}{Time} & \\multicolumn{1}{c}{Improv time} \\\\
	\t\\midrule
	\t\\endfirsthead
	\t\\caption{Results - continued}\\\\
	\t\\toprule
	\t\\multicolumn{1}{c}{\$\\begin{array}{c}
	\\text{Problem}\\\\ \\text{Size}
	\\end{array}\$} & Method & \\multicolumn{1}{c}{LDE} & \\multicolumn{1}{c}{RMSD} & \\multicolumn{1}{c}{Time} & \\multicolumn{1}{c}{Improv time} \\\\
	\t\\midrule
	\t\\endhead")

	c = 1e-9
	for prob in list_of_problems
		# before make the benchmark we need to verify if is possible solve the problem
		# considering the limit_time
		data = preprocessing(string("HCProtInsctances/",string(prob,".nmr")))
		try 
			solq = conformation(data,optq,limit_time)
			solc = conformation(data,optc,limit_time)
		catch 
			print(" 💥 The problem $(prob) cannot be solved within the timeout.\n")
			print(" \n")
			println(io, prob )
            println(ioperf,"$(prob) , Inf, Inf ")
			continue
		end
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
		LDEc = outputfilter(solc,"lde")
		PTc = f(bchc).time*c
		
		improv = (-1.0+PTc/PTq)*100
		rmsd = evalrmsd(solc,solq)[2]

        println(ioperf,"$(prob) , $(PTc), $(PTq) ")

		println(iogen,"$(prob) & Classic & $(@sprintf("%.3e",LDEc)) &  & $(@sprintf("%.4e",PTc)) & \\\\")
		println(iogen,"$(data.dim) & Quaternion & $(@sprintf("%.3e",LDEq)) & $(@sprintf("%.4e",rmsd)) & $(@sprintf("%.4e",PTq)) & $(@sprintf("%1.5f",improv))\\\\ \\cline{2-6} \\addlinespace")
	end

	println(iogen, "\\end{xltabular}")

	close(io)
    close(ioperf)
	close(iogen)
end

function performRMSD(limit_time,opwrite::String="a",f::Function=median)
	# setup 
	#optq = ConformationSetup(1.0e-8,quaternionBP,true,true,1.0e-8)
	optc = ConformationSetup(1.0e-8,classicBP,true,true,1.0e-8)
	data = preprocessing("toyinstance.nmr")
	# first run (to optimize)
	#solq = conformation(data,optq,limit_time)
	solc = conformation(data,optc,limit_time)
	#list_of_problems = ["pdb1a57", "pdb1b4c", "pdb1ba5"]
	list_of_problems = ["pdb1a57", "pdb1b4c", "pdb1ba5", "pdb1d1n", "pdb1dp3", "pdb1du1", "pdb1eii", "pdb1fcl", "pdb1fd6", "pdb1hf9", "pdb1i2u", "pdb1i2v", "pdb1ijc", "pdb1jlz", "pdb1jw3", "pdb1k0v", "pdb1k2h", "pdb1k36", "pdb1k37", "pdb1kuw", "pdb1kz0", "pdb1kz2", "pdb1kz5", "pdb1lvz", "pdb1m4e", "pdb1ma2", "pdb1ma4", "pdb1ma5", "pdb1ma6", "pdb1mpe", "pdb1nd9", "pdb1ne5", "pdb1nmj", "pdb1o53", "pdb1oqp", "pdb1plw", "pdb1plx", "pdb1pv0", "pdb1pzr", "pdb1qlk", "pdb1r57", "pdb1ry3", "pdb1s4h", "pdb1s4j", "pdb1s6j", "pdb1sa8", "pdb1t2y", "pdb1t5q", "pdb1tot", "pdb1v6r", "pdb1v92", "pdb1vd7", "pdb1vd9", "pdb1vdb", "pdb1vpc", "pdb1wnk", "pdb1wo4", "pdb1wo5", "pdb1x60", "pdb1x9v", "pdb1y5c", "pdb1yx7", "pdb1yx8", "pdb1yxr", "pdb2a2y", "pdb2a4j", "pdb2adl", "pdb2adn", "pdb2ajj", "pdb2ajm", "pdb2ajn", "pdb2ajo", "pdb2akk", "pdb2bzb", "pdb2c0s", "pdb2dci", "pdb2eem", "pdb2fva", "pdb2fvf", "pdb2fxz", "pdb2g9j", "pdb2g9l", "pdb2h5m", "pdb2hep", "pdb2j0z", "pdb2j10", "pdb2j11", "pdb2jmy", "pdb2jn5", "pdb2jnk", "pdb2jpn", "pdb2jua", "pdb2jvd", "pdb2jwe", "pdb2jws", "pdb2jwu", "pdb2jxf", "pdb2jz5", "pdb2k2a", "pdb2k2f", "pdb2k36", "pdb2k37", "pdb2k3i", "pdb2k6s", "pdb2k7o", "pdb2kbm", "pdb2kdh", "pdb2kdl", "pdb2kdm", "pdb2kdp", "pdb2kdr", "pdb2kes", "pdb2kib", "pdb2kjn", "pdb2kjo", "pdb2kjr", "pdb2kko", "pdb2kl5", "pdb2klz", "pdb2ko1", "pdb2koz", "pdb2kp0", "pdb2ksg", "pdb2kt8", "pdb2kuh", "pdb2kwh", "pdb2kxa", "pdb2kyb", "pdb2l3m", "pdb2l3n", "pdb2l45", "pdb2l5r", "pdb2l6q", "pdb2l6r", "pdb2l98", "pdb2lci", "pdb2lde", "pdb2le2", "pdb2le7", "pdb2ler", "pdb2lgi", "pdb2lhc", "pdb2lhd", "pdb2lhe", "pdb2lhg", "pdb2lix", "pdb2lld", "pdb2lm9", "pdb2lmf", "pdb2ln3", "pdb2lo2", "pdb2lqp", "pdb2lr0", "pdb2lrh", "pdb2ls9", "pdb2lsa", "pdb2lss", "pdb2lt3", "pdb2lu6", "pdb2lwa", "pdb2lx0", "pdb2lxb", "pdb2lxy", "pdb2lxz", "pdb2lzf", "pdb2m1a", "pdb2m1j", "pdb2m2y", "pdb2m3f", "pdb2m4j", "pdb2m8m", "pdb2m8o", "pdb2m97", "pdb2m9g", "pdb2m9r", "pdb2mdv", "pdb2me1", "pdb2me2", "pdb2me3", "pdb2me4", "pdb2mg1", "pdb2mg2", "pdb2mg3", "pdb2mh8", "pdb2mhw", "pdb2mi1", "pdb2mi7", "pdb2mij", "pdb2mix", "pdb2mj1", "pdb2mj2", "pdb2mjg", "pdb2mji", "pdb2mle", "pdb2mlf", "pdb2mo5", "pdb2mpu", "pdb2mpz", "pdb2mra", "pdb2msu", "pdb2mty", "pdb2muh", "pdb2mvj", "pdb2mvt", "pdb2mvx", "pdb2mz6", "pdb2n00", "pdb2n0b", "pdb2n0c", "pdb2n0d", "pdb2n0e", "pdb2n0f", "pdb2n0g", "pdb2n0h", "pdb2n35", "pdb2n3a", "pdb2n41", "pdb2n4e", "pdb2n5q", "pdb2n67", "pdb2n6n", "pdb2n7i", "pdb2n7j", "pdb2n9b", "pdb2na9", "pdb2ndc", "pdb2nde", "pdb2ndk", "pdb2nvj", "pdb2p6j", "pdb2p81", "pdb2pqe", "pdb2rlg", "pdb2rlh", "pdb2rmy", "pdb2rnd", "pdb2roo", "pdb2rql", "pdb2rut", "pdb2ruv", "pdb2rux", "pdb2ruy", "pdb2rv3", "pdb2rv5", "pdb2y4q"]

	io = open("not_solved.csv",opwrite)
	iogen = open("table_general.tex",opwrite)
	#ioop  = open("table_operations.tex",opwrite)

	println(iogen, "\\begin{xltabular}{\\textwidth}{r|rS[table-format=1.3e+2]S[table-format=1.4e+2]S[table-format=1.4e+2]S[table-format=-1.5]}
	\t\\caption{Results} \\label{tab:genResults}\\\\
	\t\\toprule
	\t\\multicolumn{1}{c}{Problem} & Method & \\multicolumn{1}{c}{LDE} & \\multicolumn{1}{c}{RMSD} & \\multicolumn{1}{c}{Time} & \\multicolumn{1}{c}{Improv time} \\\\
	\t\\midrule
	\t\\endfirsthead
	\t\\caption{Results - continued}\\\\
	\t\\toprule
	\t\\multicolumn{1}{c}{Problem} & Method & \\multicolumn{1}{c}{LDE} & \\multicolumn{1}{c}{RMSD} & \\multicolumn{1}{c}{Time} & \\multicolumn{1}{c}{Improv time} \\\\
	\t\\midrule
	\t\\endhead")

	# println(ioop, "\\begin{xltabular}{\\textwidth}{r|rcccSc}
	# \t\\caption{Operation Results}\\label{tab:opResults}\\\\
	# \t\\toprule
	# \t\\multicolumn{1}{c}{Problem} & Method & \\multicolumn{1}{c}{\\makecell[cc]{Operation \\\\Nodes}} & \\multicolumn{1}{c}{\\makecell[cc]{Operation \\\\Virtual Path}} & \\multicolumn{1}{c}{\\makecell[cc]{Operation \\\\ddf}} & \\multicolumn{1}{c}{\\makecell[cc]{No. Operations\\\\ Improvments}} & \\multicolumn{1}{c}{\\makecell[cc]{Branchs \\\\ Prunes}}  \\\\
	# \t\\midrule
	# \t\\endfirsthead
	# \t\\caption{Results - continued}\\\\
	# \t\\toprule
	# \t\\multicolumn{1}{c}{Problem} & Method & \\multicolumn{1}{c}{\\makecell[cc]{Operation \\\\Nodes}} & \\multicolumn{1}{c}{\\makecell[cc]{Operation \\\\Virtual Path}} & \\multicolumn{1}{c}{\\makecell[cc]{Operation \\\\ddf}} &
	# \t\\multicolumn{1}{c}{\\makecell[cc]{No. Operations\\\\ Improvments}} & \\multicolumn{1}{c}{\\makecell[cc]{Branchs \\\\ Prunes}}  \\\\
	# \t\\midrule
	# \t\\endhead")
	
	c = 1e-9
	for prob in list_of_problems
		# before make the benchmark we need to verify if is possible solve the problem
		# considering the limit_time
		data = preprocessing(string("HCProtInsctances/",string(prob,".nmr")))
		try 
			#solq = conformation(data,optq,limit_time)
			solc = conformation(data,optc,limit_time)
		catch 
			print(" 💥 The problem $(prob) cannot be solved within the timeout.\n")
			print(" \n")
			println(io, prob )
			continue
		end
		print(" 🔔 The problem $(prob) can be solved within the timeout. \n") 
		print(" 👏 The benchmark will be performed for both methods! \n")
		# print(" 🕐 Running benchmark using QuaternionBP... \n")
		# bchq = @benchmark conformation($(data),$(optq),$(limit_time)+Second(10.0))
		# print(" 🏁 The benchmark using QuaternionBP in $(prob) was done!\n")
		# LDEq = outputfilter(solq,"lde")
		# PTq = f(bchq).time*c
 
		print(" 🕞 Running benchmark using ClassicBP... \n")
		bchc = @benchmark conformation($(data),$(optc),$(limit_time)+Second(10.0))
		print(" 🏁 The benchmark using ClassicBP in $(prob) was done!\n")
		print(" \n")
		LDEc = outputfilter(solc,"lde")
		PTc = f(bchc).time*c
		
		#improv = (-1.0+PTc/PTq)*100

		println(iogen,"$(prob) & Classic & $(@sprintf("%.3e",LDEc)) &  & $(@sprintf("%.4e",PTc)) & \\\\")

		originalcoord = originalxyz(prob)
		coordsol = outputfilter(solc,"xyz")
		rmsdval = Array{Float64,1}(undef,size(coordsol)[2])
		for i = 1:size(coordsol)[2]
			rmsdval[i] = rmsd(vec2array(coordsol[:,i]),originalcoord)[2]
			println(iogen,"&  &  & $(@sprintf("%.4e",rmsdval[i])) & & \\\\")
		end
		
		println(iogen,"& minimum &  & $(@sprintf("%.4e",minimum(rmsdval))) & & $(argmin(rmsdval)) \\\\ \\cline{2-6} \\addlinespace")

		#rmsd = evalrmsd(solc,solq)[2]

		#println(iogen,"& Quaternion & $(@sprintf("%.3e",LDEq)) & $(@sprintf("%.4e",rmsd)) & $(@sprintf("%.4e",PTq)) & $(@sprintf("%1.5f",improv))\\\\ \\cline{2-6} \\addlinespace")

		#improv_nop = (-1.0 .+ solc.nop.node./solq.nop.node).*100
		#improv_total = (-1.0 + sum(solc.nop.node)/sum(solq.nop.node))*100
		#improv_total = (-1.0 + (solc.nop.node[1]+solc.nop.node[2]*5+solc.nop.node[3]*16+solc.nop.node[1]*21)/(solq.nop.node[1]+solq.nop.node[2]*5+solq.nop.node[3]*16+solq.nop.node[1]*21))*100
		# println(ioop,"$(prob) & Classic & $(solc.nop.node[1]) & $(solc.nop.virtual_path[1]) & $(solc.nop.ddf[1]) & & $(solc.nop.branch) \\\\")
		# println(ioop,"& & $(solc.nop.node[2]) & $(solc.nop.virtual_path[2]) & $(solc.nop.ddf[2]) & & $(solc.nop.prune) \\\\")
		# println(ioop,"& & $(solc.nop.node[3]) & $(solc.nop.virtual_path[3]) & $(solc.nop.ddf[3]) & & \\\\")
		# println(ioop,"& & $(solc.nop.node[4]) & $(solc.nop.virtual_path[4]) & $(solc.nop.ddf[4]) & & \\\\")
		
		# println(ioop,"& Quaternion & $(solq.nop.node[1]) & $(solq.nop.virtual_path[1]) & $(solq.nop.ddf[1]) & $(@sprintf("%.2f",improv_nop[1])) & \\\\")
		# println(ioop,"& & $(solq.nop.node[2]) & $(solq.nop.virtual_path[2]) & $(solq.nop.ddf[2]) & $(@sprintf("%.2f",improv_nop[2])) & \\\\")
		# println(ioop,"& & $(solq.nop.node[3]) & $(solq.nop.virtual_path[3]) & $(solq.nop.ddf[3]) & $(@sprintf("%.2f",improv_nop[3])) & \\\\")
		# println(ioop,"& & $(solq.nop.node[4]) & $(solq.nop.virtual_path[4]) & $(solq.nop.ddf[4]) & $(@sprintf("%.2f",improv_nop[4])) & \\\\ \\cline{2-7}\\addlinespace")
	end

	println(iogen, "\\end{xltabular}")
	# println(ioop, "\\end{xltabular}")

	close(io)
	close(iogen)
	# close(ioop)
end

function performoneRMSD(prob, bpall=false, limit_time=Second(120.0); ε=1.0e-8, ε_virtual=1.0e-8)
    optc = ConformationSetup(ε, classicBP, bpall, true, ε_virtual)
    data = preprocessing("toyinstance.nmr")

    global solc = conformation(data, optc, limit_time)
    
	c = 1e-9
    # before make the benchmark we need to verify if is possible solve the problem
    # considering the limit_time
    data = preprocessing(string("HCProtInsctances/", string(prob, ".nmr")))
    try
        solc = conformation(data, optc, limit_time)
    catch
        print(" 💥 The problem $(prob) cannot be solved within the timeout.\n")
        print(" \n")
		return
    end
    print(" 🔔 The problem $(prob) can be solved within the timeout. \n")
    print(" 👏 The benchmark will be performed for both methods! \n")
    
    print(" 🕞 Running benchmark using ClassicBP... \n")
    bchc = @benchmark conformation($(data), $(optc), $(limit_time) + Second(10.0))
    print(" 🏁 The benchmark using ClassicBP in $(prob) was done!\n")
    print(" \n")``
    LDEc = outputfilter(solc, "lde")
    #PTc = f(bchc).time * c

	@show LDEc
	#@show PTc

    originalcoord = originalxyz(prob)
    coordsol = outputfilter(solc, "xyz")
    rmsdval = Array{Float64,1}(undef, size(coordsol)[2])
    for i = 1:size(coordsol)[2]
        rmsdval[i] = rmsd(vec2array(coordsol[:, i]), originalcoord)[2]
    end
	@show rmsdval

	@show (argmin(rmsdval), minimum(rmsdval))
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

function originalxyz(prob)
	b = readdlm(string("HCProtInsctances/",string(prob,".xyz")))
	return Float64.(b[2:length(b[:,1]),2:4])
	#return [A[i,:] for i in 1:length(A[:,1])]
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