using MolecularConformation,BenchmarkTools, Plots, DataFrames, CSV, Printf, Dates

gr()

function plot_results(df::DataFrame;kargs...)
	gdf = groupby(df,:method)
	for info in [:mean,:median,:minimum,:maximum]
		#index = [j for j in 1:length(gdf[1].ref) if gdf[1].size[j]<=gdf[i].size[1]]
		plot(gdf[2].size[1:end],gdf[2][1:end,info],color = [:black],line = (:dot,1),xaxis= ("Size of problems"),yaxis =("Mean processing time (s)"), label = "QuaternionBP", yformatter = :scientific, legend=:topleft; kargs...);
		plot!(gdf[1].size[1:end],gdf[1][1:end,info],color = [:black],label = "ClassicBP");
		savefig("results/figures/perf_$(info)_$(df.ref[1]).pdf");
		#png("results/figures/perf_$(info)_$(df.ref[1])");
	end
end

function replot(ndiag=[3,4,5,10,100,200,300,400,500,600,700,800,900,1000]; infos::Array{Symbol,1}=[:mean,:median,:minimum,:maximum], fixedsize=false, fixedindex=false, improv=false, sample=true, improvFunction::Function = (PTq,PTc) -> (-1.0+PTc/PTq)*100, memory::Bool=false, kargs...)
	if fixedsize
		df = DataFrame()
		for diag in ndiag
			df = vcat(df, CSV.read("results/$(diag).csv",DataFrame))
		end
		gdf = groupby(df,[:size,:method])
		for i in 1:2:length(gdf)	
			for info in infos
				if fixedindex 
					index = 1:length(gdf[i].ref)
				else
					index = [j for j in 1:length(gdf[i].ref) if gdf[i].ref[j]<=gdf[i].size[1]]
				end
				if memory
					plot(gdf[i].ref[index],map(improvFunction, gdf[i+1][index,:memory], gdf[i][index,:memory]),color = [:black],xaxis= ("Amount of extra distances to each vertex"),yaxis =("Memory utilization ratio"), legend = false; kargs...);
					if sample
						plot!(twinx(),gdf[i].ref[index],gdf[i][index,:samples],color = [:green],label = "Benchmark Samples");
					end
					savefig("results/figures/perf_$(info)_size$(gdf[i].size[1])_memory.pdf");
				end
				if improv
					plot(gdf[i].ref[index],map(improvFunction, gdf[i+1][index,info], gdf[i][index,info]),color = [:black],xaxis= ("Amount of extra distances to each vertex"),yaxis =("Improvement percentage"), legend = false; kargs...);
					if sample
						plot!(twinx(),gdf[i].ref[index],gdf[i][index,:samples],color = [:green],label = "Benchmark Samples");
					end
					savefig("results/figures/perf_$(info)_size$(gdf[i].size[1])_improv.pdf");
				else
					plot(gdf[i+1].ref[index],map(x -> x/gdf[i].size[1], gdf[i+1][index,info]),color = [:black],line = (:dot,1),xaxis= ("Amount of extra distances"),yaxis =("Mean processing time (s)"), label = "QuaternionBP", yformatter = :scientific, legend=:topleft; kargs...);
					plot!(gdf[i].ref[index],map(x -> x/gdf[i].size[1], gdf[i][index,info]),color = [:black],label = "ClassicBP");
					if sample
						plot!(twinx(),gdf[i+1].ref[index],gdf[i+1][index,:samples],color = [:green],label = "Benchmark Samples Q");
						plot!(twinx(),gdf[i].ref[index],gdf[i][index,:samples],color = [:red],label = "Benchmark Samples C");
					end

					savefig("results/figures/perf_$(info)_size$(gdf[i].size[1]).pdf");
				end
			end
		end
	else
		for diag in ndiag
			if improv
				df = CSV.read("results/$(diag).csv",DataFrame)
				gdf = groupby(df,:method)
				for info in [:mean,:median,:minimum,:maximum]
					plot(gdf[2].size[1:end],map(improvFunction, gdf[2][1:end,info], gdf[1][1:end,info]),color = [:black],xaxis= ("Size of problems"),yaxis =("Improvement percentage"), legend = false; kargs...);
					savefig("results/figures/perf_$(info)_$(df.ref[1]).pdf");
				end
			else
				df = CSV.read("results/$(diag).csv",DataFrame)
				plot_results(df; kargs...)
			end
		end	
	end
end

function stringtovec(s::Union{String,Array{Int64,1}})
	if typeof(s)==Array{Int64,1}
		return s
	end
	vn = Int[]
	p = 2
	for k=2:length(s)-1
			if s[k]==','
				push!(vn,parse(Int,s[p:k-1]))
				p = k+1
			end
	end
	push!(vn,parse(Int,s[p:length(s)-1]))
	return vn
end

function write_results(df::DataFrame; improv::Function = (PTc, PTq) -> (-1.0+PTc/PTq)*100)
	gdf = groupby(df,:method)
	n = length(gdf[1].mean)
	# average percentual imṕrovement in processing time
	improv_geometric_mean = improv((prod(gdf[1].mean))^(1/n), (prod(gdf[2].mean))^(1/n))
	improv_mean = sum(improv.(gdf[1].mean, gdf[2].mean))/n
	improv_median = sum(improv.(gdf[1].median,gdf[2].median))/n
	improv_minimum = sum(improv.(gdf[1].minimum,gdf[2].minimum))/n
	improv_maximum = sum(improv.(gdf[1].maximum,gdf[2].maximum))/n
	samples_c = gdf[1].samples
	samples_q = gdf[2].samples
	io = open("results/improvs_$(df.ref[1]).txt", "w");
	write(io, "Remark: Analysis of improvements from quaternion to classic \n")
	write(io, "Average of improvements in processing time (%)\n")
	write(io, "mean -> $(improv_mean) \n");
	write(io, "median -> $(improv_median) \n");
	write(io, "minimum -> $(improv_minimum) \n");
	write(io, "maximum -> $(improv_maximum) \n");
	write(io, "geometric mean over means -> $(improv_geometric_mean) \n");
	write(io, "samples classic -> $(samples_c) \n");
	write(io, "samples quaternion -> $(samples_q) \n");
	close(io)
end

function rewrite(ndiag=[3,4,5,10,100,200,300,400,500,600,700,800,900,1000])
	for diag in ndiag
		df = CSV.read("results/$(diag).csv",DataFrame)
		write_results(df)
	end
end

function createtableimprov(ndiag=[3,4,5,10,100,200,300,400,500,600,700,800,900,1000]; improv::Function = (PTq,PTc) -> (-1.0+PTc/PTq)*100, ispercent=false)
	headerText = ""
	geometric_mean1Text = ""
	geometric_mean2Text = ""
	improvText = ""
	comma = ""
	for diag in ndiag
		gdf = groupby(CSV.read("results/$(diag).csv",DataFrame),:method)
				
		headerText = string(headerText, comma, "\\multicolumn{1}{c|}{$(diag)}")
		n = length(gdf[1].mean)
		
		geometric_mean1 = (prod(gdf[1].mean))^(1/n)
		geometric_mean1Text = string(geometric_mean1Text, comma, @sprintf("%.2e", geometric_mean1))
		
		geometric_mean2 = (prod(gdf[2].mean))^(1/n)
		geometric_mean2Text = string(geometric_mean2Text, comma, @sprintf("%.2e", geometric_mean2))
		
		improvText = string(improvText,comma, "\\multicolumn{1}{c|}{",@sprintf("%.2f", improv(geometric_mean2,geometric_mean1)),"$(ispercent ? "\\%" : "")}")
		
		comma == "" && (comma = " & ")
	end
	io = open("results/tex/improvtable.tex", "w");

	write(io, "
	\\begin{table}[H]
        \\centering
        {\\footnotesize
        
        \\begin{tabular}{||r |$(repeat("S[table-format=1.2e+2] |", length(ndiag)-1)) S[table-format=1.2e+2]||}
                \\hline
				")
	write(io, "        \$n_d\$ & $(headerText) \\\\\n");
	write(io, "        \\hline\n");
	write(io, "        Classic & $(geometric_mean1Text) \\\\\n");
	write(io, "        Quater. & $(geometric_mean2Text) \\\\\n");
	write(io, "        Improv & $(improvText) \\\\\n");
	write(io, "        \\hline
	\\end{tabular}}
	\\caption{Improvement percentage in geometric means of \\texttt{QuaternionBP} in relation to \\texttt{ClassicBP} considering results of the benchmark.}
	\\label{table:improvlavor}
\\end{table}");


	close(io)
end

"""
This function was created to make tests in order to perform our algorithms. Essentially, 
to execute this routine is necessary to type in REPL

julia> perform(d)

where d is an integer value defined by the vector [3,4,5,10,50,100,200,300,400,500,600,700,800,900,1000,2000,3000]. Warning: some files need to be download.

"""
# ndiag is defined by [3,4,5,10,50,100,200,300,400,500,600,700,800,900,1000,2000,3000]
function perform(ndiag,allsolutions=false,MDE=false; ε=1.0e-4, virtual_ε=1.0e-8,benchmarkSeconds::Int64=60, benchmarkSamples::Int64=100000, improv::Function = (PTc, PTq) -> (-1.0+PTc/PTq)*100)
	#initialization
	opt_classic = ConformationSetup(ε,classicBP,allsolutions,MDE,virtual_ε)
	opt_quaternion = ConformationSetup(ε,quaternionBP,allsolutions,MDE,virtual_ε)
	data = preprocessing("toyinstance.nmr")
	conformation(data,opt_quaternion)
	conformation(data,opt_classic)
	folders = ["10","100","200","300","400","500","600","700","800","900","1000","2000","3000"]
	c = 10.0^(-9)
	df = DataFrame(ref = Int[],method=String[],size = Int[],mean=Float64[],median=Float64[],minimum=Float64[],maximum=Float64[],samples = Int[],memory = Int[],allocs = Int[],gcmean=Float64[],gcmedian=Float64[],gcminimum=Float64[],gcmaximum=Float64[])
	print(" 🔔 Starting perform with nd = $(ndiag). Did you set up the stack limit of your OS ? I hope so!\n")
	for foldername in folders
		cd(foldername)
		psize = parse(Int,foldername)
		if ndiag > psize
			data = preprocessing(string(foldername,"-$(psize).nmr"))
		else
			data = preprocessing(string(foldername,"-$(ndiag).nmr"))
		end
		bcha = @benchmarkable conformation($(data),$(opt_classic)) seconds=benchmarkSeconds samples=benchmarkSamples
		bch = run(bcha)
		push!(df,[ndiag,"classicBP",psize,mean(bch).time*c,median(bch).time*c,minimum(bch).time*c, maximum(bch).time*c,length(bch.times),bch.memory,bch.allocs,mean(bch).gctime*c,median(bch).gctime*c,minimum(bch).gctime*c, maximum(bch).gctime*c])


		bcha = @benchmarkable conformation($(data),$(opt_quaternion)) seconds=benchmarkSeconds samples=benchmarkSamples
		bch = run(bcha)
		push!(df,[ndiag,"quaternionBP",psize,mean(bch).time*c,median(bch).time*c,minimum(bch).time*c, maximum(bch).time*c,length(bch.times),bch.memory,bch.allocs,mean(bch).gctime*c,median(bch).gctime*c,minimum(bch).gctime*c, maximum(bch).gctime*c])
		
		cd("..")
		print(" 🎉 Benchmarks using nd = $(ndiag) and the problem with $(psize) atoms were done!\n")
	#	display(df)
		write_results(df)
	end
	CSV.write("results/$(ndiag).csv", df;delim=";")
	#plot_results(df)
	write_results(df,improv=improv)
end

"""
	runperf

This function is used to run several tests

# Example
```
julia-repl
julia> runperf()

```
some graphs will be saved in current folder. 
"""
function runperf(;ndiags::Vector{Int64}=[3,4,5,10,100,200,300,400,500,600,700,800,900,1000],ε=1.0e-4, virtual_ε=1.0e-8,benchmarkSeconds::Int64=60, benchmarkSamples::Int64=100000,improv::String = "(PTc, PTq) -> (-1.0+PTc/PTq)*100")
	
	iolog = open("log.txt","w")
	println(iolog,"function_of_improvement=$(improv)")
	println(iolog,"precision=$(ε)")
	println(iolog,"precision_virtual=$(virtual_ε)")
	println(iolog,"benchmark_seconds=$(benchmarkSeconds)")
	println(iolog,"benchmark_samples_limit=$(benchmarkSamples)")
	println(iolog,"list_of_diagonals_used=$(ndiags)")
	println(iolog,"start_time=$(Dates.now())")
    close(iolog)
	
	for ndiag in ndiags # [3,10,100,400,800] 
		perform(ndiag, ε=ε, virtual_ε=virtual_ε, benchmarkSeconds=benchmarkSeconds,benchmarkSamples=benchmarkSamples,improv=eval(Meta.parse(improv)))
	end

	iolog = open("log.txt","a")
	println(iolog,"end_time=$(Dates.now())")
    close(iolog)
end
