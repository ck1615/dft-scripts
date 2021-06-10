#!/usr/bin/env python3
"""
This file contains functions for generating an input file for the correct ibrav
as the relaxed output file for PWscf calculations
"""
import sys
import getopt
import os
from ase.io import read, write
from misctools import strindex, strindices as stri
from ase.units import Bohr

def get_cmdline_options():
    """
    This function parses the command line options to get the input filename and
    any options requested
    """

    #Define PWscf calculation types: 
    pwscf_calctypes = ['scf', 'nscf', 'bands', 'relax', 'vc-relax', 'md','vc-md']

    try:
        seed = sys.argv[1]
    except IndexError:
        IndexError("No filename has been provided.")
    argv = sys.argv[2:]

    #Iterate through options
    kwargs = {}
    opts, args = getopt.getopt(argv, "c:")
    for opt, arg in opts:
        if opt in ['-c', '--calculation']:
            kwargs['calc_type'] = arg
            assert arg in pwscf_calctypes,\
            "Chosen calculation type: {}, is not one of the options for PWscf."
        else:
            Warning("Calculation type not specified, we'll assume to use the"+\
                    "same as what was in the output file.")

    return seed, kwargs

class QEOutput:


    def __init__(self, seed, calc_type='scf'):

        self.calc_type = calc_type
        self.outputfile = seed
        self.inputfile = self.outputfile.replace("out", "in")
        self.celldms = None
        self.elempos = None
        self.input_lines = None
        self.relaxed_lines = None
        self.final_fname = None
        return None


    def get_new_structure(self):

        self.Atoms = read(self.outputfile)
        self.celldms = [l/Bohr if i in range(3) else l for i,l in \
                enumerate(self.Atoms.cell.cellpar())]

        #Get elements as named in the quantum-espresso input file
        if self.input_lines is None:
            self.get_input_lines()

        old_atompos_idx = strindex(self.input_lines, "ATOMIC_POSITIONS crystal")

        #Assume ATOMIC_POSITIONS is LAST part of input file
        atompos_lines = self.input_lines[old_atompos_idx + 1:]
        elems = [line.split()[0] for line in atompos_lines]

        posns = self.Atoms.get_scaled_positions()
        self.elempos = {tuple(pos): elems[i] for i, pos in enumerate(posns)}

        return None


    def get_input_lines(self):

        with open(self.inputfile, "r+") as f:
            self.input_lines = f.readlines()

        return None


    def replace_structure(self):

        #Update lattice parameters (in Bohr)
        lattice_updated_lines = self.replace_lattice_params()

        #Update scaled positions
        relaxed_lines = self.replace_scaled_positions(lattice_updated_lines)

        return relaxed_lines


    def get_ibrav(self):
        if self.input_lines is None:
            self.get_input_lines()

        ibrav_idx = strindex(self.input_lines, "ibrav", first=True)
        self.ibrav = int(self.input_lines[ibrav_idx].split()[2].strip(","))

        return None


    def replace_lattice_params(self):

        #Get ibrav
        self.get_ibrav()

        #Get cell dimension indices required for given self.ibrav
        if self.ibrav in [6, 7] :
            idx_list = [1,3]
            self.celldms[2] /= self.celldms[0]
        elif self.ibrav == 8:
            idx_list = [1,2,3]
            self.celldms[1] /= self.celldms[0]
            self.celldms[2] /= self.celldms[0]
        elif self.ibrav == 1:
            idx_list = [1]

        for i in idx_list:
            celldm_idx = strindex(self.input_lines, "celldm({})".format(i))

            old_celldm = self.input_lines[celldm_idx].partition("celldm({})".\
                    format(i))[-1].split()[1]
            if old_celldm[-1] == ",":
                old_celldm = old_celldm[:-1]

            self.input_lines[celldm_idx] = self.input_lines[celldm_idx].\
                    replace(old_celldm, "{}d0".format(str(self.celldms[i - 1])))

        return self.input_lines


    def replace_scaled_positions(self, lattice_updated_lines):
        #Get index for atomic position
        atom_pos_idx = strindex(lattice_updated_lines, "ATOMIC_POSITIONS crystal")

        #Assume ATOMIC_POSITIONS is LAST part of input file
        lattice_updated_lines = lattice_updated_lines[:atom_pos_idx + 1]

        for pos, elem in self.elempos.items():
            lattice_updated_lines.append('{} {:.15f} {:.15f} {:.15f}\n'.format(elem, pos[0],\
                    pos[1], pos[2]))

        return lattice_updated_lines

    def modify_calc_type(self, final_lines):

        idx_calc_type = strindex(final_lines, "calculation=")
        final_lines[idx_calc_type] = "  calculation='{}'\n".format(self.calc_type)

        #Change name to reflect calculation name
        init_calc_type = self.inputfile.split(".")[-2]
        self.final_fname = self.inputfile.replace(init_calc_type,
                "relaxed.{}".format(self.calc_type))

        return final_lines


    def print_relaxed_structure(self):

        #Define input file and get line
        self.get_input_lines()

        #Get celldms & elempos
        self.get_new_structure()

        #Replace structure in input_lines
        final_lines = self.replace_structure()

        #Change calculation type
        final_lines = self.modify_calc_type(final_lines)

        #Get the Mg-O bond length
        #mgo = get_mgo_bond_length(celldms, elempos)

        #Implement changes for distance-based constraints
        #final_lines = distance_constraints(relaxed_lines, mgo)


        #Print new file
        with open(self.final_fname, "w+") as f:
            for line in final_lines:
                f.write(line)

def main():

    seed, kwargs = get_cmdline_options()

    Structure = QEOutput(seed, **kwargs)
    Structure.print_relaxed_structure()

if __name__ == "__main__":
        main()


