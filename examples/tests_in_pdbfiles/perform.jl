using MolecularConformation, PrettyTables

function perform(;allsolutions=false,highlight="PT",color = :yellow)

	opt_classic = ConformationSetup(0.000001,classicBP,allsolutions)
	opt_quaternion = ConformationSetup(0.000001,quaternionBP,allsolutions)

	# primary run to optimize
	data = preprocessing("toyinstance.nmr")
	sol = conformation(data,opt_classic);
	sol = conformation(data,opt_quaternion);
	# run to all
	list_of_problems = ["pdb1a03.nmr","pdb1a57.nmr","pdb1a7f.nmr"]
	table_header = ["problem", "method", "LDE", "PT ", "Number of solutions"]
	# defing array to storage table 
	content = Array{Any,2}(undef,2*length(list_of_problems),5)
	k=1
	for prob in list_of_problems
		data = preprocessing(prob)
		content[k,2] = "quaternionBP" 
		sol = conformation(data,opt_quaternion)
		content[k,1] = prob
		content[k,3] = outputfilter(sol,"lde")
		content[k,4] = sol.elapsedtime
		content[k,5] = sol.number
		k = k+1

	
		content[k,2] = "classicBP" 
		sol = conformation(data,opt_classic)
		content[k,1] = prob
		content[k,3] = outputfilter(sol,"lde")
		content[k,4] = sol.elapsedtime
		content[k,5] = sol.number
		k = k+1
	end
	if highlight == "PT"
		H1 = Highlighter((content,i,j)->isodd(i)&&content[i,4]==minimum(content[i:i+1,4]), Crayon(foreground = color))
		H2 = Highlighter((content,i,j)->!isodd(i)&&content[i,4]==minimum(content[i-1:i,4]), Crayon(foreground = color))
	elseif	highlight == "LDE"
		H1 = Highlighter((content,i,j)->isodd(i)&&content[i,3]==minimum(content[i:i+1,3]), Crayon(foreground = color))
		H2 = Highlighter((content,i,j)->!isodd(i)&&content[i,3]==minimum(content[i-1:i,3]), Crayon(foreground = color))
	else
		error("highlight option not recognized")
	end
	
	pretty_table(content,table_header;highlighters=(H1,H2))
end


