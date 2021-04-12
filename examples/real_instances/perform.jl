using MolecularConformation, PrettyTables, DelimitedFiles, BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000

include("rmsd.jl")

"""
The perform functions is used to run tests in order to compare algorithms.
The use is pretty simple but depends on PrettyTables.jl package. Consequently, 
the user need to install this package (installation: pkg> add PrettyTables).
To run perform script just type (inside examples/tests_in_pdbfiles folder):

julia> include(perform.jl)
julia> perform()

As output a table is printed and .tex file is write in file. Several arguments can 
be modified. By default the perform() function has the following arguments
perform(;writefile="latex",allsolutions=false,highlight="PT",color = :yellow)

As optional argument we can setup if we want all solutions:
julia> perform(allsolutions=true)

Or what aspect of the table given by output we want highlighted:
julia> perform(highlight="PT") # for best processing time
julia> perform(highlight="LDE") #for best LDE

Changing writefile="html" and output.html file is provided.

"""
function perform(f::Function=median;writefile=:latex,allsolutions=false,highlight="PT",color = :yellow)

	opt_classic = ConformationSetup(0.000001,classicBP,allsolutions)
	#opt_classicOpt = ConformationSetup(0.000001,classicBPOpt,allsolutions)
	opt_quaternion = ConformationSetup(0.000001,quaternionBP,allsolutions)

	# primary run to optimize
	data = preprocessing("toyinstance.nmr")
	solc = conformation(data,opt_classic);
	#sol = conformation(data,opt_classicOpt);
	solq = conformation(data,opt_quaternion);
	# run to all
#	list_of_problems = ["pdb1a03","pdb1a57","pdb1a7f"]
#	list_of_problems = ["pdb1a03","pdb1a57","pdb1a7f","pdb1acz","pdb2l2g","pdb2l2i","pdb2l3b","pdb2l3d","pdb2l32","pdb2l33","pdb1bct","pdb2jmy","pdb2kxa"]
#	list_of_problems = ["pdb1a03"]
	list_of_problems = ["pdb2kxa","pdb2jmy","pdb1a57"] # Prolina free list
	table_header = ["problem", "method", "LDE", "PT ", "Num. sol","rmsd","Improv"]
	table_header2 = ["problem", "method", " [ +- , / , √ ] node", " [+- , / , √ ] virtual path"," [ +- , / , √ ] ddf ", " Num. Branch ", " Num. Pru ", "Improv. nop", "Improv. tot"]
	# defing array to storage table 
	content = Array{Any,2}(undef,2*length(list_of_problems),7)
	content2 = Array{Any,2}(undef,2*length(list_of_problems),9)
	k=1
	c = 10.0^(-9)
	for prob in list_of_problems
		data = preprocessing("$(prob).nmr")
		solq = conformation(data,opt_quaternion)
		bch = @benchmark conformation($(data),$(opt_quaternion))
		content[k,1] = prob
		content[k,2] = "quaternionBP" 
		content[k,3] = outputfilter(solq,"lde")
		content[k,4] = f(bch).time*c
		content[k,5] = solq.number
		content[k,6] = " - " 
		content2[k,1] = prob
		content2[k,2] = "quaternionBP" 
		content2[k,3] = solq.nop.node
		content2[k,4] = solq.nop.virtual_path
		content2[k,5] = solq.nop.ddf
		content2[k,6] = solq.nop.branch
		content2[k,7] = solq.nop.prune
		
		k = k+1
	
	
		content[k,2] = "classicBP" 
		solc = conformation(data,opt_classic)
		bch = @benchmark conformation($(data),$(opt_classic))
		content[k,1] = prob
		content[k,3] = outputfilter(solc,"lde")
		content[k,4] = f(bch).time*c
		content[k,5] = solc.number
		content[k,6] = evalrmsd2(solc,solq)[2]
		content[k-1,7] = (-1.0+content[k,4]/content[k-1,4])*100 #tc = tq+ptq 
		content[k,7] = " - "

		content2[k,1] = prob
		content2[k,2] = "classicBP" 
		content2[k,3] = solc.nop.node
		content2[k,4] = solc.nop.virtual_path
		content2[k,5] = solc.nop.ddf
		content2[k,6] = solc.nop.branch
		content2[k,7] = solc.nop.prune
		content2[k-1,8] = (-1.0.+content2[k,3]./content2[k-1,3]).*100
		content2[k,8] = " - "
		content2[k-1,9] = (-1.0+sum(content2[k,3])/sum(content2[k-1,3]))*100
		content2[k,9] = " - "
		k = k+1
		
	#sol = conformation(data,opt_classic)
	#bch = @benchmark conformation($(data),$(opt_classic))
	#content[k,1] = prob
	#content[k,2] = "quaternionBP" 
	#content[k,3] = outputfilter(sol,"lde")
	#content[k,4] = f(bch).time*c
	#content[k,5] = sol.number
	#content[k,6] = evalrmsd(sol,"$(prob).xyz")[2]
	#content2[k,1] = prob
	#content2[k,2] = "quaternionBP" 
	#content2[k,3] = sol.nop.node
	#content2[k,4] = sol.nop.virtual_path
	#content2[k,5] = sol.nop.ddf
	#content2[k,6] = sol.nop.branch
	#content2[k,7] = sol.nop.prune
	#k = k+1
	

	end
	#if highlight == "PT"
	#	H1 = Highlighter((content,i,j)->isodd(i)&&content[i,4]==minimum(content[i:i+1,4]), Crayon(foreground = color))
	#	H2 = Highlighter((content,i,j)->!isodd(i)&&content[i,4]==minimum(content[i-1:i,4]), Crayon(foreground = color))
	#elseif	highlight == "LDE"
	#	H1 = Highlighter((content,i,j)->isodd(i)&&content[i,3]==minimum(content[i:i+1,3]), Crayon(foreground = color))
	#	H2 = Highlighter((content,i,j)->!isodd(i)&&content[i,3]==minimum(content[i-1:i,3]), Crayon(foreground = color))
	#elseif highlight == "NOP"
	#	H3 = Highlighter((content2,i,j)->isodd(i)&&content2[i,3]==minimum(content2[i:i+1,3]), Crayon(foreground = color))
	#	H4 = Highlighter((content2,i,j)->!isodd(i)&&content2[i,3]==minimum(content2[i-1:i,3]), Crayon(foreground = color))
	#else
	#	error("highlight option not recognized")
	#end
	
	#pretty_table(content,table_header;formatters=ft_printf("%5.3f",[4,7]),highlighters=(H1,H2))
	#pretty_table(content2,table_header2;formatters=ft_printf("%5.3f",[4,7]),highlighters=(H1,H2))
	if writefile==:latex
		open("output.tex","w") do f
			pretty_table(f,content,table_header,backend=:latex)
			pretty_table(f,content2,table_header2,backend=:latex)
		end
	else	
		open("output.html","w") do f
			pretty_table(f,content,table_header,backend=:html)
			pretty_table(f,content2,table_header2,backend=:html)
		end
	end
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

function evalrmsd(sols::ConformationOutput,file::String)
	coordsols = outputfilter(sols,"xyz")
	valrmsd = [1000.0,1000.0]
	xyzfile = readdlm(file)
	coordfile = Float64.(copy(xyzfile[2:xyzfile[1,1]+1,2:4]))
	#display(coordfile)
	for k=1:sols.number
		coords = vec2array(coordsols[:,k])	
		valrmsdaux = rmsd(coords,coordfile)
		if valrmsdaux[2] <= valrmsd[2]
			valrmsd = [valrmsdaux[1],valrmsdaux[2]]
		end
	end
	return valrmsd
end

function evalrmsd2(sols1::ConformationOutput,sols2::ConformationOutput)
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
