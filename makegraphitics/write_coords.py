import os
import numpy as np


class Writer(object):
    def __init__(self, sim, system_name="comment line"):
        # Takes a numpy 3xN array of atom coordinates and outputs
        # them in different formats for viewing/modelling
        self.coords = sim.coords
        self.atom_labels = sim.atom_labels
        self.molecule = sim.molecule_labels
        self.charges = sim.atom_charges

        self.bonds = sim.bonds
        self.bond_labels = sim.bond_labels
        self.angles = sim.angles
        self.angle_labels = sim.angle_labels
        self.dihedrals = sim.dihedrals
        self.dihedral_labels = sim.dihedral_labels
        self.impropers = sim.impropers
        self.improper_labels = sim.improper_labels

        self.box_dimensions = sim.box_dimensions
        self.system_name = system_name

        if hasattr(sim, "masses"):
            self.masses = sim.masses
        if hasattr(sim, "bond_coeffs"):
            self.bond_coeffs = sim.bond_coeffs
        if hasattr(sim, "angle_coeffs"):
            self.angle_coeffs = sim.angle_coeffs
        if hasattr(sim, "dihedral_coeffs"):
            self.dihedral_coeffs = sim.dihedral_coeffs
        if hasattr(sim, "improper_coeffs"):
            self.improper_coeffs = sim.improper_coeffs
        if hasattr(sim, "pair_coeffs"):
            self.pair_coeffs = sim.pair_coeffs

        type_lists = {
            "atom": "pair_coeffs",
            "bond": "bond_coeffs",
            "angle": "angle_coeffs",
            "dihedral": "dihedral_coeffs",
            "improper": "improper_coeffs",
        }
        for type_list in type_lists:
            maxlabel = max(getattr(self, type_list + "_labels"))
            coeff = type_lists[type_list]
            if hasattr(sim, coeff):
                maxcoeff = max(getattr(self, coeff))
            else:
                maxcoeff = 0
            setattr(self, "n" + type_list + "_types", max(maxcoeff, maxlabel))

    def write_xyz(self, filename="out.xyz", option="w"):
        with open(filename, option) as outfile:
            outfile.write(str(len(self.coords)) + "\n" + self.system_name + "\n")
            for i in range(len(self.coords)):
                xyz = "{:>15.8f} {:>15.8f} {:>15.8f}".format(
                    self.coords[i][0], self.coords[i][1], self.coords[i][2]
                    )
                if hasattr(self, "masses"):
                    try:
                        mass = self.masses[self.atom_labels[i]][1]
                    except KeyError:
                        mass = 100
                    if abs(mass - 12.0) < 0.5:
                        atom_label = "C   "
                    elif abs(mass - 1.0) < 0.5:
                        atom_label = "H   "
                    elif abs(mass - 14.0) < 0.5:
                        atom_label = "N   "
                    elif abs(mass - 16.0) < 0.5:
                        atom_label = "O   "
                    elif abs(mass - 22.9) < 0.5:
                        atom_label = "Na  "
                    elif abs(mass - 40.1) < 0.5:
                        atom_label = "Ca  "
                    else:
                        atom_label = "{:<5}".format(self.atom_labels[i])
                else:
                    atom_label = "{:<5}".format(self.atom_labels[i])

                outfile.write(atom_label + xyz + "\n")
            print "Coords written to " + str(filename)

    def write_pdb(self, filename="out.pdb", resname='GRA'):
        with open(filename, "w") as outfile:
            outfile.write("TITLE     {}\n".format(self.system_name))
            outfile.write("CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} {:<11}{:>4}\n"
                          .format(
                              self.box_dimensions[0, 1] - self.box_dimensions[0, 0],
                              self.box_dimensions[1, 1] - self.box_dimensions[1, 0],
                              self.box_dimensions[2, 1] - self.box_dimensions[2, 0],
                              90.0, 90.0, 90.0, "P 1", 1
                            )
                          )

            outfile.write("MODEL        1\n")
            for i in range(len(self.coords)):
                j = self.atom_labels[i]
                if hasattr(self, "masses"):
                    try:
                        mass = self.masses[j][1]
                    except KeyError:
                        mass = 100
                    if abs(mass - 12.0) < 0.5:
                        atom_label = "C"
                    elif abs(mass - 1.0) < 0.5:
                        atom_label = "H"
                    elif abs(mass - 14.0) < 0.5:
                        atom_label = "N"
                    elif abs(mass - 16.0) < 0.5:
                        atom_label = "O"
                    elif abs(mass - 22.9) < 0.5:
                        atom_label = "Na"
                    elif abs(mass - 40.1) < 0.5:
                        atom_label = "Ca"
                    else:
                        atom_label = str(j)
                else:
                    atom_label = str(j)

                outfile.write(
                    "ATOM  {:>5} {:<4} {:<3}  {:>4}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}\n"
                    .format(i + 1, atom_label, resname, 1,
                            self.coords[i][0], self.coords[i][1], self.coords[i][2], 1.0, 0.0)
                )

            print "Coords written to " + str(filename)

    def write_reaxff(self, filename="data.lammps"):
        # atom_type charge
        masses = np.unique(self.masses.values())
        nreax_types = len(masses)
        # reax_types = {mass:i+1 for i,mass in enumerate(masses)}
        reax_types = {12.011: 1, 1.008: 2, 15.999: 3}
        with open(filename, "w") as outfile:
            outfile.write(
                "# "
                + self.system_name
                + "\n"
                + str(len(self.coords))
                + " atoms \n"
                + "\n"
                + str(nreax_types)
                + " atom types \n"
                + "\n"
                + str(self.box_dimensions[0, 0])
                + "\t"
                + str(self.box_dimensions[0, 1])
                + "\t xlo xhi \n"
                + str(self.box_dimensions[1, 0])
                + "\t"
                + str(self.box_dimensions[1, 1])
                + "\t ylo yhi \n"
                + str(self.box_dimensions[2, 0])
                + "\t"
                + str(self.box_dimensions[2, 1])
                + "\t zlo zhi \n"
                + "\n"
            )
            if hasattr(self, "masses"):
                outfile.write("\n Masses \n \n")
                for mass in reax_types:
                    outfile.write(str(reax_types[mass]) + "\t" + str(mass) + "\n")

            outfile.write("\n Atoms \n \n")

            for i in range(len(self.coords)):
                atom_type = self.atom_labels[i]
                reax_type = reax_types[self.masses[atom_type][1]]
                outfile.write(
                    str(i + 1)
                    + "\t "
                    + str(reax_type)  # atom ID
                    + "\t "
                    + str(self.charges[i])  # atom type
                    + "\t "
                    + str(self.coords[i][0])  # atomcharg
                    + "\t "
                    + str(self.coords[i][1])  # x
                    + "\t "
                    + str(self.coords[i][2])  # y
                    + "\n "  # z
                )

            print "Coords written to " + filename

    def write_lammps(self, filename="data.lammps"):
        # atom_type full
        with open(filename, "w") as outfile:
            outfile.write(
                "# "
                + self.system_name
                + "\n"
                + str(len(self.coords))
                + " atoms \n"
                + str(len(self.bonds))
                + " bonds \n"
                + str(len(self.angles))
                + " angles \n"
                + str(len(self.dihedrals))
                + " dihedrals \n"
                + str(len(self.impropers))
                + " impropers \n"
                "\n"
                + str(self.natom_types)
                + " atom types \n"
                + str(self.nbond_types)
                + " bond types \n"
                + str(self.nangle_types)
                + " angle types \n"
                + str(self.ndihedral_types)
                + " dihedral types \n"
                + str(self.nimproper_types)
                + " improper types \n"
                + "\n"
                + str(self.box_dimensions[0, 0])
                + "\t"
                + str(self.box_dimensions[0, 1])
                + "\t xlo xhi \n"
                + str(self.box_dimensions[1, 0])
                + "\t"
                + str(self.box_dimensions[1, 1])
                + "\t ylo yhi \n"
                + str(self.box_dimensions[2, 0])
                + "\t"
                + str(self.box_dimensions[2, 1])
                + "\t zlo zhi \n"
                + "\n"
            )
            if hasattr(self, "masses"):
                outfile.write("\n Masses \n \n")
                for i in self.masses:
                    outfile.write(str(i) + "\t" + str(self.masses[i][1]) + "\n")

            if hasattr(self, "pair_coeffs"):
                outfile.write("\n Pair Coeffs \n \n")
                for i in self.pair_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.pair_coeffs[i][1])
                        + "\t"
                        + str(self.pair_coeffs[i][2])
                        + "\n"
                    )

            if hasattr(self, "bond_coeffs"):
                outfile.write("\n Bond Coeffs \n \n")
                for i in self.bond_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.bond_coeffs[i][1])
                        + "\t"
                        + str(self.bond_coeffs[i][2])
                        + "\n"
                    )
            if hasattr(self, "angle_coeffs"):
                outfile.write("\n Angle Coeffs \n \n")
                for i in self.angle_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.angle_coeffs[i][1])
                        + "\t"
                        + str(self.angle_coeffs[i][2])
                        + "\n"
                    )
            if hasattr(self, "dihedral_coeffs"):
                outfile.write("\n Dihedral Coeffs \n \n")
                for i in self.dihedral_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.dihedral_coeffs[i][1])
                        + "\t"
                        + str(self.dihedral_coeffs[i][2])
                        + "\t"
                        + str(self.dihedral_coeffs[i][3])
                        + "\t"
                        + str(self.dihedral_coeffs[i][4])
                        + "\n"
                    )

            if hasattr(self, "improper_coeffs"):
                outfile.write("\n Improper Coeffs \n \n")
                for i in self.improper_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.improper_coeffs[i][1])
                        + "\t"
                        + str(self.improper_coeffs[i][2])
                        + "\n"
                    )

            outfile.write("\n Atoms \n \n")

            for i in range(len(self.coords)):
                outfile.write(
                    str(i + 1)
                    + "\t "
                    + str(self.molecule[i])  # atom ID
                    + "\t "
                    + str(self.atom_labels[i])  # molecule ID
                    + "\t "
                    + str(self.charges[i])  # atom type
                    + "\t "
                    + str(self.coords[i][0])  # atomcharg
                    + "\t "
                    + str(self.coords[i][1])  # x
                    + "\t "
                    + str(self.coords[i][2])  # y
                    + "\n "  # z
                )

            if len(self.bonds):
                outfile.write("\n Bonds \n \n")
                for i in range(len(self.bonds)):
                    outfile.write(
                        str(i + 1)
                        + "\t "
                        + str(self.bond_labels[i])  # bond ID
                        + "\t "
                        + str(self.bonds[i][0])
                        + "\t "
                        + str(self.bonds[i][1])  # atom 1
                        + "\n"  # atom 2
                    )

            if len(self.angles):
                outfile.write("\n Angles \n \n")
                for i in range(len(self.angles)):
                    outfile.write(
                        str(i + 1)
                        + "\t "
                        + str(self.angle_labels[i])  # angle ID
                        + "\t "
                        + str(self.angles[i][0])
                        + "\t "
                        + str(self.angles[i][1])
                        + "\t "
                        + str(self.angles[i][2])
                        + "\n"
                    )

            if len(self.dihedrals):
                outfile.write("\n Dihedrals \n \n")
                for i in range(len(self.dihedrals)):
                    outfile.write(
                        str(i + 1)
                        + "\t"
                        + str(self.dihedral_labels[i])
                        + " \t"
                        + str(self.dihedrals[i][0])
                        + " \t"
                        + str(self.dihedrals[i][1])
                        + " \t"
                        + str(self.dihedrals[i][2])
                        + " \t"
                        + str(self.dihedrals[i][3])
                        + " \n"
                    )

            if len(self.impropers):
                outfile.write("\n Impropers \n \n")
                for i in range(len(self.impropers)):
                    outfile.write(
                        str(i + 1)
                        + "\t"
                        + str(self.improper_labels[i])
                        + "\t"
                        + str(self.impropers[i][0])
                        + " \t"
                        + str(self.impropers[i][1])
                        + " \t"
                        + str(self.impropers[i][2])
                        + " \t"
                        + str(self.impropers[i][3])
                        + " \n"
                    )

            print "Coords written to " + filename

    def write_gromacs(self, gro_filename="gromacs.gro", itp_filename="gromacs.itp", resname="GRA"):
        with open(gro_filename, "w") as outfile:
            outfile.write(self.system_name + "\n")
            outfile.write(str(len(self.coords)) + "\n")
            for i in range(len(self.coords)):
                j = self.atom_labels[i]
                if hasattr(self, "masses"):
                    try:
                        mass = self.masses[j][1]
                    except KeyError:
                        mass = 100
                    if abs(mass - 12.0) < 0.5:
                        atom_label = "C"
                    elif abs(mass - 1.0) < 0.5:
                        atom_label = "H"
                    elif abs(mass - 14.0) < 0.5:
                        atom_label = "N"
                    elif abs(mass - 16.0) < 0.5:
                        atom_label = "O"
                    elif abs(mass - 22.9) < 0.5:
                        atom_label = "Na"
                    elif abs(mass - 40.1) < 0.5:
                        atom_label = "Ca"
                    else:
                        atom_label = str(j)
                else:
                    atom_label = str(j)

                outfile.write(
                    "{:>5}{:<5}{:<5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}\n"
                    .format(1, resname, atom_label, i + 1,
                            0.1*self.coords[i][0], 0.1*self.coords[i][1], 0.1*self.coords[i][2])
                )
            outfile.write(
                    "{:>11.5f}{:>11.5f}{:>11.5f}\n".format(
                        0.1*(self.box_dimensions[0, 1] - self.box_dimensions[0, 0]),
                        0.1*(self.box_dimensions[1, 1] - self.box_dimensions[1, 0]),
                        0.1*(self.box_dimensions[2, 1] - self.box_dimensions[2, 0])
                        )
            )

        with open(itp_filename, "w") as outfile:
            outfile.write(
                "; "
                + self.system_name
                + "\n\n"
            )
            outfile.write(
                "[ atomtypes ]\n"
                + ";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb\n"
            )
            dicts = {d[2]: [k, d[1]] for k, d in self.masses.items()}
            masses = {l: {1: m, 2: k} for k, [l, m] in dicts.items()}
            for i in masses:
                outfile.write(
                    "{:<8} {:<11} {:>8.5f} {:>8.5f} {:>5} {:>13.5e} {:>13.5e}\n"
                    .format(masses[i][2], masses[i][2], 0.0, 0.0, 'A',
                            0.1*self.pair_coeffs[i][2], 4.184*self.pair_coeffs[i][1])
                )
            outfile.write("\n")
            outfile.write(
                "[ atoms ]\n"
                + ";   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type\n"
            )
            qtot = 0.0
            for i in range(len(self.coords)):
                j = self.atom_labels[i]
                if hasattr(self, "masses"):
                    try:
                        mass = self.masses[j][1]
                    except KeyError:
                        mass = 100
                    if abs(mass - 12.0) < 0.5:
                        atom_label = "C"
                    elif abs(mass - 1.0) < 0.5:
                        atom_label = "H"
                    elif abs(mass - 14.0) < 0.5:
                        atom_label = "N"
                    elif abs(mass - 16.0) < 0.5:
                        atom_label = "O"
                    elif abs(mass - 22.9) < 0.5:
                        atom_label = "Na"
                    elif abs(mass - 40.1) < 0.5:
                        atom_label = "Ca"
                    else:
                        atom_label = str(j)
                else:
                    atom_label = str(j)

                qtot += self.charges[i]
                outfile.write("{:>6} {:>4} {:>5} {:>5} {:>5} {:>4} {:>12.6f} {:>12.5f} ; qtot {:>6.3f}\n"
                    .format(i + 1, self.masses[j][2], self.molecule[i], resname, atom_label,
                            i + 1, self.charges[i], self.masses[j][1], qtot)
                )
            outfile.write("\n")
            outfile.write(
                "[ bonds ]\n"
                + ";   ai     aj funct   r             k\n"
            )
            for i in range(len(self.bonds)):
                j = self.bond_labels[i]
                iat, jat = self.bonds[i]
                iatl = self.atom_labels[iat-1]; jatl = self.atom_labels[jat-1]
                outfile.write(
                    "{:>6} {:>6} {:>3} {:>13.4e} {:>14.5e} ; {:>6} - {:>6}\n".format(
                        iat, jat, 1,
                        0.1*self.bond_coeffs[j][2], 836.8*self.bond_coeffs[j][1],
                        self.masses[iatl][2], self.masses[jatl][2]
                    )
                )
            outfile.write("\n")
            outfile.write(
                "[ angles ]\n"
                + ";   ai     aj     ak    func    theta         cth\n"
            )
            for i in range(len(self.angles)):
                j = self.angle_labels[i]
                iat, jat, kat = self.angles[i]
                iatl = self.atom_labels[iat-1]; jatl = self.atom_labels[jat-1]; katl = self.atom_labels[kat-1]
                outfile.write(
                    "{:>6} {:>6} {:>6} {:>6} {:>13.4e} {:>14.5e} ; {:>6} - {:>6} - {:>6}\n"
                    .format(iat, jat, kat, 1,
                            self.angle_coeffs[j][2], 8.368*self.angle_coeffs[j][1],
                            self.masses[iatl][2], self.masses[jatl][2], self.masses[katl][2]
                            )
                )
            outfile.write("\n")
            outfile.write(
                "[ dihedrals ] ; propers\n"
                + ";    i      j      k      l   func   c1           c2           c3           c4\n"
            )
            for i in range(len(self.dihedrals)):
                j = self.dihedral_labels[i]
                iat, jat, kat, lat = self.dihedrals[i]
                iatl = self.atom_labels[iat-1]; jatl = self.atom_labels[jat-1];
                katl = self.atom_labels[kat-1]; latl = self.atom_labels[lat-1]
                outfile.write(
                    "{:>6} {:>6} {:>6} {:>6} {:>6} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e} ; {:>6} - {:>6} - {:>6} - {:>6}\n"
                    .format(iat, jat, kat, lat, 5,
                            4.184*self.dihedral_coeffs[j][1], 4.184*self.dihedral_coeffs[j][2],
                            4.184*self.dihedral_coeffs[j][3], 4.184*self.dihedral_coeffs[j][4],
                            self.masses[iatl][2], self.masses[jatl][2],
                            self.masses[katl][2], self.masses[latl][2]
                            )
                )
            outfile.write("\n")
            outfile.write(
                "[ dihedrals ] ; impropers\n"
                + ";    i      j      k      l   func    theta         cth\n"
            )
            for i in range(len(self.impropers)):
                j = self.improper_labels[i]
                iat, jat, kat, lat = self.impropers[i]
                iatl = self.atom_labels[iat-1]; jatl = self.atom_labels[jat-1];
                katl = self.atom_labels[kat-1]; latl = self.atom_labels[lat-1]
                outfile.write(
                    "{:>6} {:>6} {:>6} {:>6} {:>6} {:>12.6e} {:>12.6e} ; {:>6} - {:>6} - {:>6} - {:>6}\n"
                    .format(iat, jat, kat, lat, 2,
                            self.improper_coeffs[j][2], 8.368*self.improper_coeffs[j][1],
                            self.masses[iatl][2], self.masses[jatl][2],
                            self.masses[katl][2], self.masses[latl][2]                           
                            )
                )

            print "Coords written to   " + gro_filename
            print "Topology written to " + itp_filename

    def write_gromacs_charmm(self, gro_filename="gromacs.gro", itp_filename="gromacs.itp", resname="GRA"):

        """
        # all OPLS atom types that are introduced by oxidation
        sim.vdw_defs = {
            1: 90,  # Cg, graphitic (aromatic)
            2: 91,  # Hg, graphitic edge
            3: 101,  # Ct, tertiary C-OH (Alcohol R3COH)
            4: 96,  # Oa, C-OH (Alcohol -OH)
            5: 97,  # Ha, C-OH (Alcohol -OH)
            6: 122,  # Oe, epoxy (Dialkyl Ether -O-)
            11: 108,  # Cb, Benzyl (Phenol C-OH)
            7: 109,  # Oa, C-OH (Phenol -OH)
            8: 209,  # Cc, Carboxylic carbon (Carboxylic Acid -COOH)
            9: 210,  # Oc, Ketone oxygen (Carboxylic Acid C=O)
            10: 211,  # Oa, alcohol (Carboxylic Acid -OH)
            12: 213, # C, carboxylate -COO (Carboxylate COO-)
            13: 214, # O, carboxylate -COO (Carboxylate COO-)
            14: 101, # Ct, tertiary C-OH (Alcohol R3COH)
            349: 349, # Na+
            354: 354, # Ca 2+ 
        }  # OPLS definitions
        masses = {1:   {1: 12.011, 2: 'CA'},
                  2:   {1: 1.008,  2: 'HA'},
                  3:   {1: 12.011, 2: 'CT'},
                  4:   {1: 15.999, 2: 'OH'},
                  5:   {1: 1.008,  2: 'HO'},
                  6:   {1: 15.999, 2: 'OS'},
                  7:   {1: 15.999, 2: 'OH'},
                  8:   {1: 12.011, 2: 'C'},
                  9:   {1: 15.999, 2: 'O'},
                  10:  {1: 15.999, 2: 'OH'},
                  11:  {1: 12.011, 2: 'CA'},
                  12:  {1: 12.011, 2: 'C'},
                  13:  {1: 15.999, 2: 'O2'}
                  14:  {1: 12.011, 2: 'CT'},
                  349: {1: 22.99,  2: 'Na'},
                  354: {1: 40.08,  2: 'Ca'}}
        """
        # CGenFF definitions (parameter from DOI: 10.1016/j.commatsci.2015.12.030)
        masses = {1:   {1: 12.011, 2: 'CGOA'},
                  2:   {1: 1.008,  2: 'HGOA'},
                  3:   {1: 12.011, 2: 'CGOOH'},
                  4:   {1: 15.999, 2: 'OGOOH'},
                  5:   {1: 1.008,  2: 'HGOOH'},
                  6:   {1: 15.999, 2: 'OGOEPO'},
                  7:   {1: 15.999, 2: 'OGOPOH'},
                  8:   {1: 12.011, 2: 'CGOCA'},
                  9:   {1: 15.999, 2: 'OGOCA2'},
                  10:  {1: 15.999, 2: 'OGOCA1'},
                  11:  {1: 12.011, 2: 'CGOPOH'},
                  12:  {1: 12.011, 2: 'CGOCB'},
                  13:  {1: 15.999, 2: 'OGOCB'},
                  14:  {1: 12.011, 2: 'CGOEPO'},
                  349: {1: 22.99,  2: 'SOD'},
                  354: {1: 40.08,  2: 'CAL'}}
        
        with open(gro_filename, "w") as outfile:
            outfile.write(self.system_name + "\n")
            outfile.write(str(len(self.coords)) + "\n")
            for i in range(len(self.coords)):
                j = self.atom_labels[i]
                if hasattr(self, "masses"):
                    try:
                        mass = self.masses[j][1]
                    except KeyError:
                        mass = 100
                    if abs(mass - 12.0) < 0.5:
                        atom_label = "C"
                    elif abs(mass - 1.0) < 0.5:
                        atom_label = "H"
                    elif abs(mass - 14.0) < 0.5:
                        atom_label = "N"
                    elif abs(mass - 16.0) < 0.5:
                        atom_label = "O"
                    elif abs(mass - 22.9) < 0.5:
                        atom_label = "Na"
                    elif abs(mass - 40.1) < 0.5:
                        atom_label = "Ca"
                    else:
                        atom_label = str(j)
                else:
                    atom_label = str(j)

                outfile.write(
                    "{:>5}{:<5}{:<5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}\n"
                    .format(1, resname, atom_label, i + 1,
                            0.1*self.coords[i][0], 0.1*self.coords[i][1], 0.1*self.coords[i][2])
                )
            outfile.write(
                    "{:>11.5f}{:>11.5f}{:>11.5f}\n".format(
                        0.1*(self.box_dimensions[0, 1] - self.box_dimensions[0, 0]),
                        0.1*(self.box_dimensions[1, 1] - self.box_dimensions[1, 0]),
                        0.1*(self.box_dimensions[2, 1] - self.box_dimensions[2, 0])
                        )
            )

        with open(itp_filename, "w") as outfile:
            outfile.write(
                "; "
                + self.system_name
                + "\n\n"
            )
            outfile.write(
                "[ moleculetype ]\n"
                + "; name  nrexcl\n"
                + resname + "     1\n"
            )
            outfile.write("\n")
            outfile.write(
                "[ atoms ]\n"
                + ";   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type\n"
            )
            qtot = 0.0
            for i in range(len(self.coords)):
                j = self.atom_labels[i]
                if hasattr(self, "masses"):
                    try:
                        mass = self.masses[j][1]
                    except KeyError:
                        mass = 100
                    if abs(mass - 12.0) < 0.5:
                        atom_label = "C"
                    elif abs(mass - 1.0) < 0.5:
                        atom_label = "H"
                    elif abs(mass - 14.0) < 0.5:
                        atom_label = "N"
                    elif abs(mass - 16.0) < 0.5:
                        atom_label = "O"
                    elif abs(mass - 22.9) < 0.5:
                        atom_label = "Na"
                    elif abs(mass - 40.1) < 0.5:
                        atom_label = "Ca"
                    else:
                        atom_label = str(j)
                else:
                    atom_label = str(j)

                qtot += self.charges[i]
                outfile.write("{:>6} {:>6} {:>5} {:>5} {:>5} {:>4} {:>12.6f} {:>12.5f} ; qtot {:>6.3f}\n"
                    .format(i + 1, masses[j][2], self.molecule[i], resname, atom_label,
                            i + 1, self.charges[i], self.masses[j][1], qtot)
                )
                #outfile.write("{:>6} {:>6} {:>5} {:>5} {:>5} {:>4} {:>12.6f} {:>12.5f} ; qtot {:>6.3f}\n"
                #    .format(i + 1, self.masses[j][2], self.molecule[i], resname, atom_label,
                #            i + 1, self.charges[i], self.masses[j][1], qtot)
                #)
            outfile.write("\n")
            outfile.write(
                "[ bonds ]\n"
                + ";   ai     aj funct   r             k\n"
            )
            for i in range(len(self.bonds)):
                iat, jat = self.bonds[i]
                iatl = self.atom_labels[iat-1]; jatl = self.atom_labels[jat-1]
                outfile.write(
                    "{:>6} {:>6} {:>3} ; {:>6} - {:>6}\n"
                    .format(iat, jat, 1, masses[iatl][2], masses[jatl][2])
                )
            outfile.write("\n")
            outfile.write(
                "[ angles ]\n"
                + ";   ai     aj     ak    func    theta         cth\n"
            )
            for i in range(len(self.angles)):
                iat, jat, kat = self.angles[i]
                iatl = self.atom_labels[iat-1]; jatl = self.atom_labels[jat-1]; katl = self.atom_labels[kat-1]
                outfile.write(
                    "{:>6} {:>6} {:>6} {:>6} ; {:>6} - {:>6} - {:>6}\n"
                    .format(iat, jat, kat, 5, masses[iatl][2], masses[jatl][2], masses[katl][2])
                )
            outfile.write("\n")
            outfile.write(
                "[ dihedrals ] ; propers\n"
                + "; ai    aj      ak      al      funct   phi0    cp      mult\n"
            )
            for i in range(len(self.dihedrals)):
                iat, jat, kat, lat = self.dihedrals[i]
                iatl = self.atom_labels[iat-1]; jatl = self.atom_labels[jat-1];
                katl = self.atom_labels[kat-1]; latl = self.atom_labels[lat-1]
                outfile.write(
                    "{:>6} {:>6} {:>6} {:>6} {:>6} ; {:>6} - {:>6} - {:>6} - {:>6}\n"
                    .format(self.dihedrals[i][0], self.dihedrals[i][1],
                            self.dihedrals[i][2], self.dihedrals[i][3], 9,
                            masses[iatl][2], masses[jatl][2], masses[katl][2], masses[latl][2])
                )
            #outfile.write("\n")
            #outfile.write(
            #    "[ dihedrals ] ; impropers\n"
            #    + "; ai    aj      ak      al      funct   phi0    cp      mult\n"
            #)
            #for i in range(len(self.impropers)):
            #    iat, jat, kat, lat = self.impropers[i]
            #    iatl = self.atom_labels[iat-1]; jatl = self.atom_labels[jat-1];
            #    katl = self.atom_labels[kat-1]; latl = self.atom_labels[lat-1]
            #    outfile.write(
            #        "{:>6} {:>6} {:>6} {:>6} {:>6} ; {:>6} - {:>6} - {:>6} - {:>6}\n"
            #        .format(self.impropers[i][0], self.impropers[i][1],
            #                self.impropers[i][2], self.impropers[i][3], 2,
            #                masses[iatl][2], masses[jatl][2], masses[katl][2], masses[latl][2])
            #    )

            print "Coords written to   " + gro_filename
            print "Topology written to " + itp_filename
