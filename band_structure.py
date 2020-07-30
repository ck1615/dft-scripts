#!/usr/bin/env python3
"""
This module extracts the values of all electronic band structure quantities
computed in a CASTEP Spectral calculation.
"""
import numpy as np
import argparse
import collections
import matplotlib.pyplot as plt
from ase.units import Hartree
import strindices as stri

def parse_seed():
    """
    This function parses the seedname for a CASTEP calculation

    Returns:
    --------
    seed: str
        The seedname of all CASTEP files
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("seed", type=str)
    seed = parser.parse_args().seed

    return seed

class BandStructure:
    """
    The Band_structure class represents all quantities describing the
    electronic band structure from a CASTEP calculation on a system.
    """

    def __init__(self):
        self.code = 'castep'
        self.nBands1 = 0            #The number of bands in spin channel 1
        self.nBands2 = 0            #The number of bands in spin channel 2 
        self.nElectrons = 0         #The number of electrons
        self.nSpins = 0             #The number of spin channels
        self.eFermi = 0             #The Fermi energy
        self.bandGap = 0            #The band gap
        self.cell = np.array([])    #The lattice parameters
        self.bgType = ""            #Direct (D) or Indirect (I) band gap
        self.eBands1 = np.array([]) #Spin 1 channel energy eigenvalues
        self.eBands2 = np.array([]) #Spin 2 channel energy eigenvalues
        self.kPoints = np.array([]) #Array of k-points (#k-points, 3)
        self.Path = np.array([])
        self.Lines = None
        self.HO = None
        self.LU = None

    def get_eFermi(self):
        for line in self.Lines:
            if self.code == 'qe':
                if 'the Fermi energy is' in line:
                    self.eFermi = float(line.split()[-2])
        return None

    def get_data(self, fname):
        with open(fname, 'r') as file:
            self.Lines = file.readlines()

    def read_bands(self, fname):
        """
        The read_bands function reads the file containing bands and extracts 
        the energy eigenvalues at each k-point, the number of bands, number of
        electrons, and the Fermi energy.

        Note for quantum-espresso verbosity must be set to 'high'.

        Parameters:
        -----------
        fname: str
            CASTEP: The name of the .bands file (<seed>.bands)
            elk: EIGVAL.OUT is assumed
            quantum-espresso: The name of the output file.
        """
        #Get data from file
        if self.Lines is None:
            self.get_data(fname)
            data = self.Lines

        #Get cell 
        self.get_cell()

        #Get bands data according to code
        if self.code == 'castep':
            self.read_castep_bands(fname, data)
        elif self.code == 'elk':
            self.read_elk_bands(data)
        elif self.code == 'qe':
            self.read_qe_bands(fname, data)
        else:
            return ValueError('Currently only CASTEP, elk or' + \
                    'quantum-espresso codes are implemented.')
        return None

    def get_cell(self):
        self.cell = np.array([[1/5.32349155533, 0.00000000000, 0.00000000000],\
            [0.00000000000, 1/5.45857015431, 0.00000000000],\
            [0.00000000000, 0.00000000000, 1/13.0877903134]])
        return None

    def get_path(self, letters=True):

        if letters:
            self.Path = {"$\Gamma$" : np.array([0.,0.,0.]), "R": \
                np.array([0.5, 0.5, 0.5]), "S": np.array([0.5, 0.5, 0]), 'T':\
                np.array([0., 0.5, 0.5]), 'U': np.array([0.5, 0., .5]), \
                "X": np.array([0.5, 0.0, 0.0]), 'Y': np.array([.0, .5, .0]),\
                "Z": np.array([0., 0., .5])}
#            self.Path = {"$\Gamma$": np.array([0.,0.,0.]), "M": np.array([0.5,\
#                    0.5, 0.5]), "$\overline{\Gamma}$": np.array([0., 1., 0.]),\
#                    "$\Gamma$'": np.array([0.,0.,0.])}
        else:
            self.Path = {"(0,0,0)" : np.array([0.,0.,0.]), "(1/2,1/2,1/2)": \
                np.array([0.5, 0.5, 0.5]), "(1/2,1/2,0)": \
                np.array([0.5, 0.5, 0]), '(0,1/2,1/2)': np.array([0., 0.5, 0.5\
                ]), '(1/2,0,1/2)': np.array([0.5, 0., .5]), "(1/2,0,0)": \
                np.array([0.5, 0.0, 0.0]), '(0,1/2,0)': np.array([.0, .5, .0]),\
                "(0,0,1/2)": np.array([0., 0., .5])}

        return None

    def read_castep_bands(self, fname, data):
        """
        This function reads a <seed>.bands file 
        from a CASTEP calculation
        """

        #Initialise dictionaries & temporary array
        e_bands_dict1 = {}
        e_bands_dict2 = {}
        k_points_dict = {}
        eK = []

        #Extracting relevant values & bands
        for iLine, line in enumerate(data):
            lineSplit = line.split()

            #Values
            if line.find("Number of spin components") > -1:
                self.nSpins = int(lineSplit[4])

            if line.find("Number of electrons") > -1:
                if self.nSpins == 1:
                    self.nElectrons = float(lineSplit[3])
                else:
                    self.nElectrons = float(lineSplit[3]) + float(lineSplit[4])

            if line.find("Number of eigenvalues") > -1:
                self.nBands1 = int(lineSplit[3])
                if self.nSpins == 2:
                    self.nBands2 = int(lineSplit[4])

            if lineSplit[0] == "Fermi":
                self.eFermi = float(lineSplit[5])*Hartree

            #K-points
            if line.find("K-point") > -1:
                k_key = int(lineSplit[1])
                k_points_dict[k_key] = np.array([float(lineSplit[2]),
                    float(lineSplit[3]), float(lineSplit[4])])

            if line.find("Spin component") > -1:
                s_key = int(lineSplit[2])

            if len(lineSplit)==1:
                eK.append(float(lineSplit[0]))
                if s_key == 1:
                    if len(eK) == self.nBands1:
                        e_bands_dict1[k_key] = eK      #Construct e_bands_dict
                        eK = []
                elif s_key == 2:
                    if len(eK) == self.nBands2:
                        e_bands_dict2[k_key] = eK      #Construct e_bands_dict
                        eK = []
                else:
                    print("Error: the number of spin components is not 1 or 2")

        #Order energy eigenvalue arrays by k-point
        e_bands_dict1 = collections.OrderedDict(sorted(e_bands_dict1.items()))
        e_bands1 = np.array(list(e_bands_dict1.values()))
        e_bands1 = e_bands1*Hartree                 #Convert eigenvalues to eV
        self.eBands1 = e_bands1

        e_bands_dict2 = collections.OrderedDict(sorted(e_bands_dict2.items()))
        e_bands2 = np.array(list(e_bands_dict2.values()))
        e_bands2 = e_bands2*Hartree          #Convert eigenvalues to eV
        self.eBands2 = e_bands2

        #Order k-point coordinates by k-point
        k_point_dict = collections.OrderedDict(sorted(k_points_dict.items()))
        self.kPoints = k_points_dict

    def read_elk_bands(self, fname):
        pass

    def read_qe_bands(self, fname, lines):
        """
        This function extracts the eigenvalues at each k-point after a nscf
        calculation performed by quantum-espresso.
        Note: verbosity of quantum-espresso calculation must be set to 'high'
        """
        #Initialise dictionaries & temporary array
        e_bands_dicts = [{},{}]
        k_points_dict = {}

        #Get indices of up- and down-spin k-points in file
        try:
            idx_up = stri.strindex(lines, '------ SPIN UP ------------')
        except:
            AssertionError("Verbosity of nscf calculation in" + \
                    " quantum-espresso must be set to 'high'")

        idx_down = stri.strindex(lines, '------ SPIN DOWN ----------')

        #Define up and down eigenvalue regions
        eigvals_up, eigvals_down = (lines[idx_up + 3:idx_down], \
                lines[idx_down + 3:])

        #Get number of electrons
        idx_nElectrons = stri.strindex(lines, 'number of electrons')
        self.nElectrons = float(lines[idx_nElectrons].split()[-1])

        #Assume spin-polarised
        self.nSpins = 2

        #Find location of k-points
        for sp, evals in enumerate([eigvals_up, eigvals_down]):
            #Define boundaries for each k-point
            kLocs = stri.strindices(evals, 'k =')
            occLocs = stri.strindices(evals, 'occupation numbers')

            #Check both numbers of k-point and occupation num. locations are 
            #equal
            assert len(kLocs) == len(occLocs), "The number of k-point and " + \
                    "'occupation number' occurences are not equal."

            #K-point eigenvalues
            for i, l in enumerate(kLocs):
                eK = []
                #Get k-point
                kps = []
                kp_string = eigvals_up[kLocs[i]].split("k =")[1].split("(")[0]
                for kval in kp_string.split():
                    degeneracy = stri.strindices(kval, '-')
                    if len(degeneracy) == 0:
                        kps.append(float(kval))
                    elif len(degeneracy) == 1:
                        if degeneracy[0] == 0:
                            kps.append(float(kval))
                        else:
                            kps.append(float(kval[:degeneracy[0]]))
                            kps.append(float(kval[degeneracy[0]:]))
                    elif len(degeneracy) == 2:
                        kps.append(float(kval[:degeneracy[1]]))
                        kps.append(float(kval[degeneracy[1]:]))
                    elif len(degeneracy) == 3:
                        kps.append(float(kval[:degeneracy[1]]))
                        kps.append(float(kval[degeneracy[1]:degeneracy[2]]))
                        kps.append(float(kval[degeneracy[2]:]))

                k_points_dict[i+1] = kps
                #Initialise
                kp_eigvals = evals[kLocs[i] + 2:occLocs[i] - 1]

                for line in kp_eigvals:
                    for val in line.split():
                        degeneracy = stri.strindices(val, '-')
                        if len(degeneracy) <= 1:
                            eK.append(float(val))
                        elif len(degeneracy) == 2:
                            eK.append(float(val[:degeneracy[1]]))
                            eK.append(float(val[degeneracy[1]:]))
                        else:
                            raise AssertionError("Too many values stuck together")

                e_bands_dicts[sp][i+1] = eK
                eK = []

            e_bands_dicts[sp] = collections.OrderedDict(sorted(e_bands_dicts[sp].items()))
            e_bands = np.array(list(e_bands_dicts[sp].values()))
            if sp == 0:
                self.eBands1 = e_bands
            elif sp ==1:
                self.eBands2 = e_bands
            else:
                raise ValueError("Spin should be 1 or 2.")

        #Order k-point coordinates by k-point
        k_point_dict = collections.OrderedDict(sorted(k_points_dict.items()))
        k_points = np.array(list(k_point_dict.values()))
        self.kPoints = k_points

    def band_gap(self):
        """
        This function gets the band numbers of the highest occupied (HO) and
        lowest unoccupied (LO) bands, computes the band gap by taking the
        difference of the maximum of the HO energies and the
        minimum of the LU energies over all k-points, and determines whether
        the band gap is direct (D) or indirect (I).
        """
        #HO & LU location
        n_lu = int(self.nElectrons/2)
        n_ho = n_lu - 1

        #Band gap & band gap type
        #Non spin-polarised case
        if self.nSpins == 1:
            self.LU = self.eBands1[:, n_lu].min()
            self.HO = self.eBands1[:, n_ho].max()
        #Spin-polarised case
        else:
            self.LU = np.array([self.eBands1[:, n_lu].min(), self.eBands2[:,\
                    n_lu].min()]).min()
            self.HO = np.array([self.eBands1[:,n_ho].max(), self.eBands2[:, \
                    n_ho].max()]).max()
 

        #Band gap
        self.bandGap = self.LU - self.HO
        #Direct band gap if indices match and indirect if not
        if self.eBands1[:, n_ho].argmax() == self.eBands1[:, n_lu].argmin():
            self.bgType = "D"
        else:
            self.bgType = "I"

        #Setting band gap to 0.00 if negative
        if self.bandGap < 0.00:
            self.bandGap = 0.00

    def plot_bs(self, fname):
        """
        This function plots the band structure using the bands parsed by the
        read_bands function.
        """
        s = 5

        #Figure title
        plotname = fname.replace(".bands", ".pdf")

        #Plot band structure
        plt.figure(1)
        plt.ylabel("$E - E_{Fermi}$ / eV")
        plt.plot(self.eBands1 - self.eFermi, "-k", linewidth=0.8)
        plt.plot(self.eBands2 - self.eFermi, "--r", linewidth=0.3)
        plt.xlim((0,len(self.eBands1 - 2)))
        plt.ylim((-7,7))
        plt.savefig(plotname)

        return

    def plot_bs_special(self, fname):

        s = 5
        #Figure title
        plotname = fname.replace(".bands", ".pdf")
        kpoints = self.kPoints
        kpoints_values = [list(item) for item in list(kpoints.values())]
        special_points = self.Path
        special_points_values = [list(item) for item in list(special_points.values())]
        special_points_keys = list(special_points.keys())

        # find indices where keys are the same 
        special_indices = [i for i, item in enumerate(kpoints_values) if item in special_points_values]

        #Plot band structure
        fig, ax = plt.subplots()

        plt.plot(self.eBands1 - self.eFermi, "-k", linewidth=0.85)
        plt.plot(self.eBands2 - self.eFermi, "--r", linewidth=0.3)

        # labels
        ax.set_ylabel(r"$E - E_{\mathrm{Fermi}}$ / eV", fontsize=12)
        ax.set_xlim((0,special_indices[-1]))
        ax.set_ylim((-8,3))

        # Plot special k points 
        ax.set_xticks(special_indices)
        ax.set_xticklabels(special_points_keys, rotation=45, fontsize=11)

        # Plot vertical lines at these points 
        for x in special_indices:
            ax.axvline(x, color="k", linewidth=0.8, alpha=0.8)

        # Plot bandgap
        #ax.axhline(self.LU - self.eFermi, color="k", linestyle="dashed", linewidth=0.8, alpha=0.8)
        #ax.axhline(self.HO - self.eFermi, color="k", linestyle="dashed", linewidth=0.8, alpha=0.8)
        plt.tight_layout()
        plt.savefig(plotname)

        return


if __name__ == "__main__":

    #Get the parsed seedname
    seed = parse_seed()
    fname = seed

    #Get the band gap
    bs = BandStructure()
    bs.read_bands(fname)
    bs.get_eFermi()
    bs.band_gap()
    bs.get_path(letters=False)
    bs.plot_bs_special(fname)

    print("Band gap: ", round(bs.bandGap, 4), " eV")
