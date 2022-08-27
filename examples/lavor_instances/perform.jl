using MolecularConformation,BenchmarkTools, Plots

gr()

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000


"""
This function was created to make tests in order to perform our algorithms. Essentially, 
to execute this routine is necessary to type in REPL

julia> perform(d)

where d is an integer value defined by the vector [3,4,5,10,50,100,200,300,400,500,600,700,800,900,1000,2000,3000]. Warning: some files need to be download.

"""
# ndiag is defined by [3,4,5,10,50,100,200,300,400,500,600,700,800,900,1000,2000,3000]
function perform(ndiag,allsolutions=false,MDE=false)
	#initialization
	opt_classic = ConformationSetup(0.000001,classicBP,allsolutions,MDE)
	opt_quaternion = ConformationSetup(0.000001,quaternionBP,allsolutions,MDE)
	data = preprocessing("toyinstance.nmr")
	sol = conformation(data,opt_quaternion)
	sol = conformation(data,opt_classic)
#	sol = conformation(data,opt_classicOpt)
	#
	folders = ["10","100","200","300","400","500","600","700","800","900","1000","2000","3000"]
	prob = parse.(Int,folders)
	tquat = zeros(length(folders),4)
	tclass = zeros(length(folders),4)
	nop  = Vector{String}(undef,length(folders))
#	tclassopt = zeros(length(folders),4)
	k = 1
	c = 10.0^(-9)
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
		tclass[k,1] = minimum(bch).time*c
		tclass[k,2] = median(bch).time*c
		tclass[k,3] = mean(bch).time*c
		tclass[k,4] = maximum(bch).time*c
	
#		bch  = @benchmark conformation($(data),$(opt_classicOpt))
#		tclassopt[k,1] = minimum(bch).time*c
#		tclassopt[k,2] = median(bch).time*c
#		tclassopt[k,3] = mean(bch).time*c
#		tclassopt[k,4] = maximum(bch).time*c


		bch = @benchmark conformation($(data),$(opt_quaternion))
		solq = conformation(data,opt_quaternion)
		nop[k] = "In $(foldername) gain was $((-1.0+sum(solc.nop.node)/sum(solq.nop.node))*100) %"
		tquat[k,1] = minimum(bch).time*c
		tquat[k,2] = median(bch).time*c
		tquat[k,3] = mean(bch).time*c
		tquat[k,4] = maximum(bch).time*c
		
		k += 1
		cd("..")
	end
	#####
	improv(q,c) = (c/q -1.0)*100
	# average improvement in minimal time
	improv_min = sum(improv.(tquat[:,1],tclass[:,1]))/length(tquat[:,1])
	# average improvement in median time
	improv_med = sum(improv.(tquat[:,2],tclass[:,2]))/length(tquat[:,2])
	# average improvement in mean time
	improv_mean = sum(improv.(tquat[:,3],tclass[:,3]))/length(tquat[:,3])
	# average improvement in maximal time
	improv_max = sum(improv.(tquat[:,4],tclass[:,4]))/length(tquat[:,4])
	io = open("improvs_$(ndiag).txt", "w");
	write(io, "minimum -> $(improv_min) \n");
	write(io, "median -> $(improv_med) \n");
	write(io, "mean -> $(improv_mean) \n");
	write(io, "maximum -> $(improv_max) \n");
	for k = 1: length(folders)
		write(io," $(nop[k]) \n")
	end
	close(io);
	#####
	j=1
	for information in ["minimum","median","mean","maximum"]
		plot(prob,tquat[:,j],color = [:black],line = (:dot,2),xaxis = ("Size of problems"),yaxis =("Average of $(information) processing time (s)"), label = "QuaternionBP")
		plot!(prob,tclass[:,j],color = [:black],label =  "ClassicBP")
		png("perf_$(information)_$(ndiag)");
		j+=1
	end	 
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
function runperf()
	for ndiag in [3,4,5,10,100,200,300,400,500,600,700,800,900,1000] # [3,10,100,200,400,800,1000] 
		perform(ndiag)
	end
end
