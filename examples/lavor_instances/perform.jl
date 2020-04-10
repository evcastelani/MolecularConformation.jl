using MolecularConformation,PyPlot,BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000


"""
This function was created to make tests in order to perform our algorithms. Essentially, 
to execute this routine is necessary to type in REPL

julia> perform(d)

where d is an integer value defined by the vector [3,4,5,10,50,100,200,300,400,500,600,700,800,900,1000,2000,3000]

"""
# ndiag is defined by [3,4,5,10,50,100,200,300,400,500,600,700,800,900,1000,2000,3000]
function perform(ndiag,allsolutions=false,LDE=false)
	#initialization
	opt_classic = ConformationSetup(0.000001,classicBP,allsolutions,LDE)
	opt_quaternion = ConformationSetup(0.000001,quaternionBP,allsolutions,LDE)
	data = preprocessing("toyinstance.nmr")
	sol = conformation(data,opt_quaternion)
	sol = conformation(data,opt_classic)
	#
	folders = ["10","100","200","300","400","500","600","700","800","900","1000","2000","3000"]
	prob = parse.(Int,folders)
	tquat = zeros(length(folders),4)
	tclass = zeros(length(folders),4)
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
		tclass[k,1] = minimum(bch).time*c
		tclass[k,2] = median(bch).time*c
		tclass[k,3] = mean(bch).time*c
		tclass[k,4] = maximum(bch).time*c

		bch = @benchmark conformation($(data),$(opt_quaternion))
		tquat[k,1] = minimum(bch).time*c
		tquat[k,2] = median(bch).time*c
		tquat[k,3] = mean(bch).time*c
		tquat[k,4] = maximum(bch).time*c
		
		k += 1
		cd("..")
	end
	j=1
	for information in ["minimum","median","mean","maximum"]
		figure()
	    	xlabel("size of problems")
		ylabel("Average of $(information) processing time (s)")
    	#grid("on")
    	#title("")
		plot(prob,tquat[:,j],"--k",label="QuaternionBP");
		plot(prob,tclass[:,j],"-k",label="ClassicBP");
	    	legend(loc="upper right",fancybox="false");
		savefig("perf_$(information)_$(ndiag).png");
		j+=1
	end	 
end
