using DelimitedFiles,MolecularConformation

opt = ConformationSetup(0.001,4.5,classical_bp,false)
qopt = ConformationSetup(0.001,4.5,quaternion_bp,false)


problems = ["inst_10.txt","inst_15.txt","inst_20.txt", "inst_30.txt","inst_40.txt","inst_50.txt","inst_60.txt","inst_70.txt","inst_80.txt",
"inst_90.txt","inst_100.txt","inst_200.txt","inst_300.txt"]

f = open("performance.tex","w")
println(f, "Method & Problem & Dim & cutoff & LDE & Time & number of solutions \\\\")

D = readdlm("inst_10.txt")
conformation(D,opt);
conformation(D,qopt);

for name in problems
    D = readdlm(name)
    (m,n) = size(D)

    s_classic = conformation(D,opt);
    s_classic = conformation(D,opt);
    println(f,"Classical_bp & $(name) & $(n) &$(opt.cutoff) & $(s_classic.molecules[1].lde) & $(s_classic.elapsedtime) & $(s_classic.number)\\\\")

    s_quaternion = conformation(D,qopt);
    s_quaternion = conformation(D,qopt);
    println(f,"Quaternion_bp& $(name) & $(n) &$(qopt.cutoff) & $(s_quaternion.molecules[1].lde) & $(s_quaternion.elapsedtime) & $(s_quaternion.number)\\\\")
end

opt = ConformationSetup(0.001,5.5,classical_bp,false)
qopt = ConformationSetup(0.001,5.5,quaternion_bp,false)


problems = ["inst_400.txt","inst_500.txt","inst_600.txt", "inst_700.txt","inst_800.txt","inst_900.txt","inst_1000.txt"]



D = readdlm("inst_10.txt")
conformation(D,opt);
conformation(D,qopt);

for name in problems
    D = readdlm(name)
    (m,n) = size(D)

    s_classic = conformation(D,opt);
    s_classic = conformation(D,opt);
    println(f,"Classical_bp & $(name) & $(n) &$(opt.cutoff) & $(s_classic.molecules[1].lde) & $(s_classic.elapsedtime) & $(s_classic.number)\\\\")

    s_quaternion = conformation(D,qopt);
    s_quaternion = conformation(D,qopt);
    println(f,"Quaternion_bp& $(name) & $(n) &$(qopt.cutoff) & $(s_quaternion.molecules[1].lde) & $(s_quaternion.elapsedtime) & $(s_quaternion.number)\\\\")
end

close(f)