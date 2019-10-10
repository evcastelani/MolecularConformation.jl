# This routine was made by Leandro Martinez https://github.com/leandromartinez98


using DelimitedFiles
include("./align.jl")

#a = readdlm("./1.xyz")
#b = readdlm("./2.xyz")

function rmsd(x,y)
	# Align x to y (same size, correspondence given by sequence)

	z = align(x,y)

	# Compute rmsd before and after alignment

	rmsd_before = 0.
	rmsd_after = 0.
	n = size(x,1)
	for i in size(x,1)
#	  global rmsd_before, rmsd_after, n
	  rmsd_before = rmsd_before + (x[i,1] - y[i,1])^2 + (x[i,2] - y[i,2])^2 + (x[i,3] - y[i,3])^2
	  rmsd_after = rmsd_after + (z[i,1] - y[i,1])^2 + (z[i,2] - y[i,2])^2 + (z[i,3] - y[i,3])^2
	end
	rmsd_before = sqrt(rmsd_before/n)
	rmsd_after = sqrt(rmsd_after/n)
	println(" rmsd before = ", rmsd_before)
	println(" rmsd after = ", rmsd_after)
	return rmsd_before,rmsd_after
end
