using MolecularConformation, PrettyTables, DelimitedFiles
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
function perform(;writefile=:latex,allsolutions=false,highlight="PT",color = :yellow)

	opt_classic = ConformationSetup(0.000001,classicBP,allsolutions)
	opt_quaternion = ConformationSetup(0.000001,quaternionBP,allsolutions)

	# primary run to optimize
	data = preprocessing("toyinstance.nmr")
	sol = conformation(data,opt_classic);
	sol = conformation(data,opt_quaternion);
	# run to all
#	list_of_problems = ["pdb1a03","pdb1a57","pdb1a7f","pdb1acz","pdb2l2g","pdb2l2i","pdb2l3b","pdb2l3d","pdb2l32","pdb2l33"]
	list_of_problems = ["pdb1a03"]
	table_header = ["problem", "method", "LDE", "PT ", "Num. sol","Num. Op.","rmsd"]
	# defing array to storage table 
	content = Array{Any,2}(undef,2*length(list_of_problems),7)
	k=1
	for prob in list_of_problems
		data = preprocessing("$(prob).nmr")
		content[k,2] = "quaternionBP" 
		sol = conformation(data,opt_quaternion)
		content[k,1] = prob
		content[k,3] = outputfilter(sol,"lde")
		content[k,4] = sol.elapsedtime
		content[k,5] = sol.number
		content[k,6] = sol.nop
		content[k,7] = evalrmsd(sol,"$(prob).xyz")
		k = k+1

	
		content[k,2] = "classicBP" 
		sol = conformation(data,opt_classic)
		content[k,1] = prob
		content[k,3] = outputfilter(sol,"lde")
		content[k,4] = sol.elapsedtime
		content[k,5] = sol.number
		content[k,6] = sol.nop
		content[k,7] = evalrmsd(sol,"$(prob).xyz")
		k = k+1
	end
	if highlight == "PT"
		H1 = Highlighter((content,i,j)->isodd(i)&&content[i,4]==minimum(content[i:i+1,4]), Crayon(foreground = color))
		H2 = Highlighter((content,i,j)->!isodd(i)&&content[i,4]==minimum(content[i-1:i,4]), Crayon(foreground = color))
	elseif	highlight == "LDE"
		H1 = Highlighter((content,i,j)->isodd(i)&&content[i,3]==minimum(content[i:i+1,3]), Crayon(foreground = color))
		H2 = Highlighter((content,i,j)->!isodd(i)&&content[i,3]==minimum(content[i-1:i,3]), Crayon(foreground = color))
	elseif highlight == "NOP"
		H1 = Highlighter((content,i,j)->isodd(i)&&content[i,6]==minimum(content[i:i+1,6]), Crayon(foreground = color))
		H2 = Highlighter((content,i,j)->!isodd(i)&&content[i,6]==minimum(content[i-1:i,6]), Crayon(foreground = color))
	else
		error("highlight option not recognized")
	end
	
	pretty_table(content,table_header;formatter=ft_printf("%5.3f",[4,7]),highlighters=(H1,H2))
	if writefile==:latex
		open("output.tex","w") do f
			pretty_table(f,content,table_header,backend=:latex)
		end
	else	
		open("output.html","w") do f
			pretty_table(f,content,table_header,backend=:html)
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
	coordfile = readdlm(file)
	for k=1:sols.number
		coords = vec2array(coordsols[:,k])	
		valrmsdaux = rmsd(coords,coordfile)
		if valrmsdaux[2] <= valrmsd[2]
			valrmsd = [valrmsdaux[1],valrmsdaux[2]]
		end
	end
	return valrmsd
end
