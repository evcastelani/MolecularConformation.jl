using DelimitedFiles,MolecularConformation, PyPlot

include("rmsd.jl")

function tests(diag)
opt = ConformationSetup(0.001,4.5,classical_bp,false)
qopt = ConformationSetup(0.001,4.5,quaternion_bp,false)


problems = ["inst_10.txt","inst_15.txt","inst_20.txt", "inst_30.txt","inst_40.txt","inst_50.txt","inst_60.txt","inst_70.txt","inst_80.txt",
"inst_90.txt","inst_100.txt","inst_200.txt","inst_300.txt"]

f = open("performance_$(diag).tex","w")
println(f, "Method & Problem & Dim & cutoff & LDE & Time & number of solutions \\\\")

prob = zeros(length(problems))
tquat = zeros(length(problems))
tclass = zeros(length(problems))

D = readdlm("inst_10.txt")
conformation(D,opt,ndiag=diag);
conformation(D,qopt,ndiag=diag);

k=1
for name in problems
    D = readdlm(name)
    (m,n) = size(D)
    prob[k] = n

    s_quaternion = conformation(D,qopt,ndiag=diag);
    println(f,"Quaternion_bp& $(name) & $(n) &$(qopt.cutoff) & $(s_quaternion.molecules[1].lde) & $(s_quaternion.elapsedtime) & $(s_quaternion.number)\\\\")
    tquat[k] = s_quaternion.elapsedtime 

    s_classic = conformation(D,opt,ndiag=diag);
    println(f,"Classical_bp & $(name) & $(n) &$(opt.cutoff) & $(s_classic.molecules[1].lde) & $(s_classic.elapsedtime) & $(s_classic.number)\\\\")
    tclass[k] = s_classic.elapsedtime

    if tquat[k]<tclass[k]
	    println(f,"Quaternion_bp was $(100*(1-tquat[k]/tclass[k])) percent faster than Classical_bp")
    else
	    println(f,"Classical_bp was $(100*(1-tclass[k]/tquat[k])) percent faster than Quaternion_bp")
    end

    k = k+1
end

opt = ConformationSetup(0.001,5.5,classical_bp,false)
qopt = ConformationSetup(0.001,5.5,quaternion_bp,false)


problems = ["inst_400.txt","inst_500.txt","inst_600.txt", "inst_700.txt","inst_800.txt","inst_900.txt","inst_1000.txt","inst_1100.txt","inst_1200.txt","inst_1300.txt","inst_1400.txt","inst_1500.txt","inst_2000.txt","inst_3000.txt","inst_4000.txt","inst_5000.txt","inst_6000.txt"]



prob = append!(prob,zeros(length(problems)))
tquat = append!(tquat,zeros(length(problems)))
tclass = append!(tclass,zeros(length(problems)))

D = readdlm("inst_10.txt")
conformation(D,opt,ndiag=diag);
conformation(D,qopt,ndiag=diag);

for name in problems
    D = readdlm(name)
    (m,n) = size(D)

    prob[k] = n
    
    s_quaternion = conformation(D,qopt,ndiag=diag);
    println(f,"Quaternion_bp& $(name) & $(n) &$(qopt.cutoff) & $(s_quaternion.molecules[1].lde) & $(s_quaternion.elapsedtime) & $(s_quaternion.number)\\\\")
    tquat[k] = s_quaternion.elapsedtime 

    s_classic = conformation(D,opt,ndiag=diag);
    println(f,"Classical_bp & $(name) & $(n) &$(opt.cutoff) & $(s_classic.molecules[1].lde) & $(s_classic.elapsedtime) & $(s_classic.number)\\\\")
    tclass[k] = s_classic.elapsedtime

    if tquat[k]<tclass[k]
	    println(f,"Quaternion_bp was $(100*(1-tquat[k]/tclass[k])) percent faster than Classical_bp")
    else
	    println(f,"Classical_bp was $(100*(1-tclass[k]/tquat[k])) percent faster than Quaternion_bp")
    end
    k = k+1
end
println(f,"sum of elapsed time for all problems in Classical_bp = $(sum(tclass))")
println(f,"sum of elapsed time for all problems in Quaternion_bp = $(sum(tquat))")
close(f)

#opt = ConformationSetup(0.001,5.5,classical_bp,false)
#qopt = ConformationSetup(0.001,5.5,quaternion_bp,false)

#s_quaternion = conformation(D,qopt);
#s_classic = conformation(D,opt);

xlabel("Size of problems ")
ylabel("Average processing time (s)")
#grid("on")
#title("")
plot(prob,tquat,"--k",label="QuaternionBP")
plot(prob,tclass,"-k",label="ClassicBP")
legend(loc="upper right",fancybox="false")
savefig("perfomance_$(diag).png")
end

function rmsd_run()

	opt = ConformationSetup(0.001,5.5,classical_bp,false)
	qopt = ConformationSetup(0.001,5.5,quaternion_bp,false)
	
	problems = ["pdb1a03v1","pdb1a03v2","pdb1a57v1","pdb1a7fv1","pdb1a7fv2","pdb1aczv1","pdb1aczv2","pdb2l2iv1","pdb2l3bv1","pdb2l3dv1","pdb2l2gv1","pdb2l2gv2","pdb2l32v1","pdb2l32v2","pdb2l33v1","pdb2l33v2"]
	
	# unsolved pdb1a57v2,pdb2l2iv2,pdb2l3bv2,pdb2l3dv2
	#
	f = open("rmsd_output.tex","w")	
	println(f,"problem  &  method & rmsd_before & rmsd_after")
	for name in problems
		print("\n Solving the problem $(name) ...\n")
		D = readdlm("PDB/$(name).txt")
		y = readdlm("PDB/$(name).xyz")
		s_classic = conformation(D,opt);
		n = length(s_classic.molecules[1].atoms)
		s = zeros(n,3)
		for i=1:n
			s[i,1] = s_classic.molecules[1].atoms[i].x
			s[i,2] = s_classic.molecules[1].atoms[i].y
			s[i,3] = s_classic.molecules[1].atoms[i].z
		end
		rb,ra = rmsd(s,y)  	
		println(f, "$(name) & BP Classic & $(rb) & $(ra)")
		s = zeros(n,3)
		s_quaternion = conformation(D,qopt);
		for i=1:n
			s[i,1] = s_quaternion.molecules[1].atoms[i].x
			s[i,2] = s_quaternion.molecules[1].atoms[i].y
			s[i,3] = s_quaternion.molecules[1].atoms[i].z
		end
		rb,ra = rmsd(s,y)
		println(f, "$(name) & BP Quaternion & $(rb) & $(ra)")
		
	end
	close(f)

end
