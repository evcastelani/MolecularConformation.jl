using MolecularConformation,BenchmarkTools, Plots, DataFrames, CSV

gr()

BenchmarkTools.DEFAULT_PARAMETERS.samples = 5

function plot_results(df::DataFrame)
	gdf = groupby(df,:method)
	for info in [:mean,:median,:minimum,:maximum]
		plot(gdf[2].size[1:end-2],gdf[2][1:end-2,info],color = [:black],line = (:dot,1),xaxis= ("Size of problems"),yaxis =("Average of mean processing time (s)"), label = "QuaternionBP", yformatter = :scientific);
		plot!(gdf[1].size[1:end-2],gdf[1][1:end-2,info],color = [:black],label = "ClassicBP");
		savefig("results/figures/perf_$(info)_$(df.ref[1]).pdf");
		#png("results/figures/perf_$(info)_$(df.ref[1])");
	end
end

function replot(ndiag=[3,4,5,10,100,200,300,400,500,600,700,800,900,1000])
	for diag in ndiag
		df = CSV.read("results/$(diag).csv",DataFrame)
		plot_results(df)
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

function write_results(df::DataFrame)
	improv(q,c) = (c/q -1.0)*100
	gdf = groupby(df,:method)
	n = length(gdf[1].mean)
	# average percentual imṕrovement in processing time
	improv_mean = sum(improv.(gdf[2].mean,gdf[1].mean))/n
	improv_median = sum(improv.(gdf[2].median,gdf[1].median))/n
	improv_minimum = sum(improv.(gdf[2].minimum,gdf[1].minimum))/n
	improv_maximum = sum(improv.(gdf[2].maximum,gdf[1].maximum))/n
	io = open("results/improvs_$(df.ref[1]).txt", "w");
	write(io, "Remark: Analysis of improvements from quaternion to classic \n")
	write(io, "Average of improvements in processing time (%)\n")
	write(io, "mean -> $(improv_mean) \n");
	write(io, "median -> $(improv_median) \n");
	write(io, "minimum -> $(improv_minimum) \n");
	write(io, "maximum -> $(improv_maximum) \n");
	close(io)
end


"""
This function was created to make tests in order to perform our algorithms. Essentially, 
to execute this routine is necessary to type in REPL

julia> perform(d)

where d is an integer value defined by the vector [3,4,5,10,50,100,200,300,400,500,600,700,800,900,1000,2000,3000]. Warning: some files need to be download.

"""
# ndiag is defined by [3,4,5,10,50,100,200,300,400,500,600,700,800,900,1000,2000,3000]
function perform(ndiag,allsolutions=false,MDE=false; ε=1.0e-4, virtual_ε=1.0e-8)
	#initialization
	opt_classic = ConformationSetup(ε,classicBP,allsolutions,MDE,virtual_ε)
	opt_quaternion = ConformationSetup(ε,quaternionBP,allsolutions,MDE,virtual_ε)
	data = preprocessing("toyinstance.nmr")
	solq = conformation(data,opt_quaternion)
	solc = conformation(data,opt_classic)
	folders = ["10","100","200","300","400","500","600","700","800","900","1000","2000","3000"]
	c = 10.0^(-9)
	df = DataFrame(ref = Int[],method=String[],size = Int[],mean=Float64[],median=Float64[],minimum=Float64[],maximum=Float64[])
	print(" 🔔 Starting perform with nd = $(ndiag). Did you set up the stack limit of your OS ? I hope so!\n")
	for foldername in folders
		cd(foldername)
		psize = parse(Int,foldername)
		if ndiag > psize
			data = preprocessing(string(foldername,"-$(psize).nmr"))
		else
			data = preprocessing(string(foldername,"-$(ndiag).nmr"))
		end
		bch  = @benchmark conformation($(data),$(opt_classic))
		solc =  conformation(data,opt_classic)
		push!(df,[ndiag,"classicBP",psize,mean(bch).time*c,median(bch).time*c,minimum(bch).time*c, maximum(bch).time*c])

		bch = @benchmark conformation($(data),$(opt_quaternion))
		solq = conformation(data,opt_quaternion)
		push!(df,[ndiag,"quaternionBP",psize,mean(bch).time*c,median(bch).time*c,minimum(bch).time*c, maximum(bch).time*c])
		
		cd("..")
		print(" 🎉 Benchmarks using nd = $(ndiag) and the problem with $(psize) atoms were done!\n")
	#	display(df)
		write_results(df)
	end
	CSV.write("results/$(ndiag).csv", df;delim=";")
	plot_results(df)
	write_results(df)
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
function runperf(;ε=1.0e-4, virtual_ε=1.0e-8)
	for ndiag in [3,4,5,10,100,200,300,400,500,600,700,800,900,1000] # [3,10,100,400,800] 
		perform(ndiag, ε=1.0e-4, virtual_ε=1.0e-8)
	end
end
