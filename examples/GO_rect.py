import makegraphitics as mg

motif = mg.molecules.Rectangle_Graphene(50, 50)
flake = mg.Crystal(motif, [1, 1, 1])

oxidiser = mg.reactors.Oxidiser(
    ratio=2.5, video_xyz=20, new_island_freq=1e14, method="rf"
)
flake = oxidiser.react(flake)

mg.Parameterise(flake, flake.vdw_defs)

name = "GO_rest"
output = mg.Writer(flake, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
output.write_reaxff(name + "reax.data")
output.write_pdb(name + ".pdb", "GO")
output.write_gromacs(name + ".gro", name + ".itp", "GO")
