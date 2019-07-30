using DelimitedFiles
@testset "Testing main functions" begin


D = readdlm("inst_10.txt")


options = ConformationSetup(0.001,4.5,classical_bp,true)

sol = conformation(D,options)


@test sol.number == 2

mol1 = MoleculeType([AtomType(0.0, 0.0, 0.0), AtomType(-1.526, 0.0, 0.0), AtomType(-2.03539, 1.43847, 0.0), AtomType(-1.46692, 2.18016, -1.2064), AtomType(-1.93283, 1.49818, -2.48956), AtomType(-1.44723, 2.29654, -3.696), AtomType(0.07718, 2.36468, -3.68164), AtomType(0.543635, 3.06946, -2.41106), AtomType(0.0165284, 2.32061, -1.19038), AtomType(0.482983, 3.0254, 0.0802008)], 4.812186979389608e-29)


@test sol.molecules[1]≈mol1


mol2 = MoleculeType(AtomType[AtomType(0.0, 0.0, 0.0), AtomType(-1.526, 0.0, 0.0), AtomType(-2.03539, 1.43847, 0.0), AtomType(-1.46692, 2.18016, 1.2064), AtomType(-1.93283, 1.49818, 2.48956), AtomType(-1.44723, 2.29654, 3.696), AtomType(0.07718, 2.36468, 3.68164), AtomType(0.543635, 3.06946, 2.41106), AtomType(0.0165284, 2.32061, 1.19038), AtomType(0.482983, 3.0254, -0.0802008)], 4.812186979389608e-29)

@test sol.molecules[2]≈mol2

options = ConformationSetup(0.001,4.5,quaternion_bp,true)

sol = conformation(D,options)


@test sol.number == 2

mol1 = MoleculeType([AtomType(0.0, 0.0, 0.0), AtomType(-1.526, 0.0, 0.0), AtomType(-2.03539, 1.43847, 0.0), AtomType(-1.46692, 2.18016, -1.2064), AtomType(-1.93283, 1.49818, -2.48956), AtomType(-1.44723, 2.29654, -3.696), AtomType(0.07718, 2.36468, -3.68164), AtomType(0.543635, 3.06946, -2.41106), AtomType(0.0165284, 2.32061, -1.19038), AtomType(0.482983, 3.0254, 0.0802008)], 4.812186979389608e-29)


@test sol.molecules[1]≈mol1


mol2 = MoleculeType(AtomType[AtomType(0.0, 0.0, 0.0), AtomType(-1.526, 0.0, 0.0), AtomType(-2.03539, 1.43847, 0.0), AtomType(-1.46692, 2.18016, 1.2064), AtomType(-1.93283, 1.49818, 2.48956), AtomType(-1.44723, 2.29654, 3.696), AtomType(0.07718, 2.36468, 3.68164), AtomType(0.543635, 3.06946, 2.41106), AtomType(0.0165284, 2.32061, 1.19038), AtomType(0.482983, 3.0254, -0.0802008)], 4.812186979389608e-29)

@test sol.molecules[2]≈mol2



end
