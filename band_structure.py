#!/usr/bin/env python3
"""
This module extracts the values of all electronic band structure quantities
computed in a CASTEP Spectral calculation.

Options to parse:
    -p --path: symbols for BZ band structure path
    -bg --band-gap: print band gap (true, false)

"""
import numpy as np
from numpy.linalg import norm
import argparse
import collections
import matplotlib.pyplot as plt
from ase.units import Hartree, Bohr
from misctools import strindices as stri
from misctools import strindex

plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=14
plt.rcParams['font.size']=16

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

def get_cmdline_options():
    """
    This function parses the command line options to establish the program and
    the type of plot
    """

    import sys
    import getopt

    #Collect filename arguments
    try:
        seed = sys.argv[1]
    except IndexError:
        IndexError("No filename has been provided.")
    argv = sys.argv[2:]

    #Iterate through options
    kwargs = {}
    opts, args = getopt.getopt(argv, "p:bg:s:c:")
    for opt, arg in opts:
        if opt in ['-p', '--path']:
            kwargs['path'] = arg
        elif opt in ['-bg', '--band-gap']:
            printBG = arg.lower()
            if printBG == "true":
                printBG = True
            else:
                printBG == False
            kwargs['printBG'] = printBG
        elif opt in ['-s', '--spin-polarised']:
            if arg.lower() == "true":
                spin_polarised = True
            else:
                spin_polarised = False
            kwargs['spinPol'] = spin_polarised
        elif opt in ['-c', '--code']:
            if arg.lower() == "castep":
                kwargs['code'] = 'castep'
            elif arg.lower() == 'q-e':
                kwargs['code'] = 'q-e'
            else:
                ValueError("Specified program {} is not castep or q-e"
                        .format(arg))

    return seed, kwargs

class BandStructure:
    """
    The Band_structure class represents all quantities describing the
    electronic band structure from a CASTEP calculation on a system.
    """


    def __init__(self, fname, printBG=True,
            code='q-e', spinPol=True):

        self.fname = fname
        self.code = code
        self.spinPolarised = spinPol
        self.calcType = None
        self.nBands1 = 0            #The number of bands in spin channel 1
        self.nBands2 = 0            #The number of bands in spin channel 2 
        self.nElectrons = 0         #The number of electrons
        self.eFermi = 0             #The Fermi energy
        self.bandGap = 0            #The band gap
        self.cell = np.array([])    #The lattice parameters
        self.bgType = ""            #Direct (D) or Indirect (I) band gap
        self.eBands1 = np.array([]) #Spin 1 channel energy eigenvalues
        self.eBands2 = np.array([]) #Spin 2 channel energy eigenvalues
        self.kPoints = {} #Array of k-points (#k-points, 3)
        self.Path = np.array([])
        self.Lines = None
        self.HO = None
        self.LU = None
        self.printBG = printBG
        self.path = None
        #Number of spin channels
        if self.spinPolarised:
            self.nSpins = 2
        else:
            self.nSpins = 1

        #Get data
        self.get_data()

        #Get alat
        self.get_alat()

        #Get cell
        self.get_cell()

        #Read bands
        self.read_bands()
        self.get_kpath_indices()
        self.get_path()
        self.get_eFermi()

    def get_alat(self):
        """
        This function gets the value of alat
        """
        self.alat = float(self.Lines[strindex(self.Lines, \
                "lattice parameter (alat)")].split()[-2]) * Bohr
        return None


    def get_eFermi(self):
        """
        This function gathers the Fermi energy in eV.
        """
        if self.code == 'q-e':
            #Open the *.scf.out file
            scf_name = self.fname.replace("bands", "scf")
            try:
                with open(scf_name, "r+") as f:
                    scfLines = f.readlines()
            except FileNotFoundError:
                with open(self.fname.replace("bands", "vc-relax"), 'r+') as f:
                    scfLines = f.readlines()

            #Assume fixed occupations and using highest occupied
            self.eFermi = float(scfLines[strindex(scfLines, 'highest occupied')].split()\
                    [-2])

        else:
            print("Fermi Energy extraction not necessary or not implemented.")
            pass
        return None


    def get_data(self):

        with open(self.fname, 'r') as file:
            self.Lines = file.readlines()


    def get_qe_calc_type(self):

        self.calcType = self.fname.split(".")[-2]


    def read_bands(self):
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
            self.get_data()

        #Get bands data according to code
        if self.code == 'castep':
            self.seed = self.fname.replace(".bands" ,"")
            self.read_castep_bands(self.Lines)
        elif self.code == 'elk':
            self.read_elk_bands(self.Lines)
        elif self.code == 'q-e':
            #Get calculation type
            self.get_qe_calc_type()
            self.read_qe_bands()
        else:
            return ValueError('Currently only CASTEP, elk or' + \
                    'quantum-espresso codes are implemented.')
        return None


    def get_cell(self):

        if self.code == 'q-e':
            self.cell = np.zeros((3,3))
            for i in range(3):
                self.cell[i] = [float(x)*self.alat for x in self.Lines[strindex(\
                    self.Lines, "a({})".format(i+1))].split()[3:6]]
        else:
            pass


    def get_path(self):

        if self.code == 'q-e':
            self.input_name = self.fname.replace("out", "in")
            with open(self.input_name, 'r+') as f:
                self.input_lines = f.readlines()

            start_idx = strindex(self.input_lines, "K_POINTS crystal_b")
            end_idx = strindex(self.input_lines, "ATOMIC_POSITIONS")
            self.path = np.array([[float(l) for l in line.split()[:-1]] for line \
                    in self.input_lines[start_idx+2:end_idx] if line.split() \
                    !=[]])


        return None


    def read_castep_bands(self, data):

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
        self.eBands1 = e_bands1 - self.eFermi       #Shift Fermi Level to 0 eV 

        e_bands_dict2 = collections.OrderedDict(sorted(e_bands_dict2.items()))
        e_bands2 = np.array(list(e_bands_dict2.values()))
        e_bands2 = e_bands2*Hartree                 #Convert eigenvalues to eV
        self.eBands2 = e_bands2  - self.eFermi      #Shift Fermi level to 0 eV

        #Order k-point coordinates by k-point
        k_point_dict = collections.OrderedDict(sorted(k_points_dict.items()))
        self.kPoints = k_points_dict

        return None


    def read_elk_bands(self, fname):
        pass


    def get_number_electrons(self):

        if self.code == "q-e":
            idx_nElectrons = strindex(self.Lines, 'number of electrons')
            self.nElectrons = float(self.Lines[idx_nElectrons].split()[4][:-1])


    def read_qe_bands(self):
        """
        This function extracts the eigenvalues at each k-point after a nscf
        calculation performed by quantum-espresso.
        Note: verbosity of quantum-espresso calculation must be set to 'high'
        """
        if self.calcType == "vc-relax":
            start_idx = strindex(self.Lines, "Begin final coordinates")
            self.Lines = self.Lines[start_idx:]
        elif self.calcType == "bands":
            start_idx = strindex(self.Lines, "End of band structure calculation")

        #Get the Fermi energy
        self.get_eFermi()

        #Get number of electrons
        self.get_number_electrons()

        #Treat spin polarised and unpolarised eigenvalue cases
        if self.nSpins == 2:
            self.get_qe_spinpol_bands()
        elif nSpins == 1:
            self.get_qe_nonspinpol_bands()

        ##Re-order eigenvalues and k-points and shift by Fermi Energy
        #e_bands_dicts[sp] = collections.OrderedDict(sorted(e_bands_dicts[sp].items()))
        #e_bands = np.array(list(e_bands_dicts[sp].values()))

        #if sp == 0:
        #    self.eBands1 = e_bands - self.eFermi
        #elif sp == 1:
        #    self.eBands2 = e_bands - self.eFermi
        #else:
        #    raise ValueError("Spin should be 1 or 2.")

        ##Order k-point coordinates by k-point
        #k_point_dict = collections.OrderedDict(sorted(k_points_dict.items()))
        #k_points = np.array(list(k_point_dict.values()))
        #self.kPoints = k_points

        return None


    def get_qe_spinpol_bands(self):
        """
        Function returns the KS eigenvalues at all k-points at the end of a
        scf, nscf, bands, or vc-relax (scf) calculation.

        Returns:
        --------
            e_bands_dicts: list of two dictionaries
                Dictionary mapping k-point to list of eigenvalues for spin-up
                and spin-down channels

            k_points_dict: dictionary mapping index to k-point vector in
            fractional coordinates
        """

        #Initialise dictionaries & temporary array
        e_bands_dicts = [{},{}]
        k_points_dict = {}

        #Get indices of up- and down-spin k-points in file
        try:
            idx_up = strindex(self.Lines, '------ SPIN UP ------------')
            idx_down = strindex(self.Lines, '------ SPIN DOWN ----------')
        except:
            AssertionError("Verbosity of nscf or bands calculation in" + \
                    " quantum-espresso must be set to 'high'")

        #Define up and down eigenvalue regions
        try:
            idx_end = strindex(self.Lines, 'the spin up/dw Fermi energies')
        except:
            try:
                idx_end = strindex(self.Lines, 'highest occupied, lowest unoccupied level (ev)')
            except:
                idx_end = strindex(self.Lines,'Writing output data file')

        eigvals_up, eigvals_down = (self.Lines[idx_up + 3:idx_down], \
                self.Lines[idx_down + 3:idx_end - 1])

        #Find location of k-points
        for sp, evals in enumerate([eigvals_up, eigvals_down]):
            kLocs, occLocs = self.get_kpoint_locations(evals)

            #K-point eigenvalues
            for i, l in enumerate(kLocs):
                #Get k-point
                k_points_dict[i+1] = self.get_kpoint(evals, l)

                #Get eigenvalues
                if self.calcType in ['scf', 'nscf', 'vc-relax']:
                    kp_eigvals = evals[kLocs[i] + 2:occLocs[i] - 1]
                elif self.calcType == 'bands':
                    if i+1 == len(kLocs):
                        kp_eigvals = evals[kLocs[i] + 2:idx_end - 1]
                    else:
                        kp_eigvals = evals[kLocs[i] + 2:kLocs[i+1]  - 1]

                e_bands_dicts[sp][i+1] = self.get_eigenvalues(kp_eigvals)

            #Re-order eigenvalues and k-points and shift by Fermi Energy
            e_bands_dicts[sp] = collections.OrderedDict(sorted(e_bands_dicts[sp].items()))
            e_bands = np.array(list(e_bands_dicts[sp].values()))

            if sp == 0:
                self.eBands1 = e_bands - self.eFermi
            elif sp == 1:
                self.eBands2 = e_bands - self.eFermi
            else:
                raise ValueError("Spin should be 1 or 2.")

            #Order k-point coordinates by k-point
            k_point_dict = collections.OrderedDict(sorted(k_points_dict.items()))
            k_points = np.array(list(k_point_dict.values()))
            self.kPoints['cartesian'] = k_points * 2*np.pi / self.alat
            self.kPoints['crystal'] = (np.matmul(k_points,\
                    self.cell) / self.alat).round(2)

        return None


    def get_kpath_indices(self):
        """
        This function takes a set of k-points and associates to it a list of
        non-negative real numbers corresponding to the distance from the
        initial k-point on the path.
        """
        path = []
        for i, kpt in enumerate(self.kPoints['cartesian']):
            if i == 0:
                path.append(0.0)
            else:
                path.append(path[i-1] + norm(kpt - \
                        self.kPoints['cartesian'][i-1]))

        #Normalise list between 0.0 and 1.0
        self.kpath_idx = [(idx - path[0]) / (path[-1] - \
                path[0]) for idx in path]

        return None


    def get_kpoint_locations(self, evals):

        """
        Returns locations of k-points for a given section of the output file.

        Optionally returns 'occupation number' location.
        """

        #Define boundaries for each k-point
        kLocs = stri(evals, 'k =')

        if self.calcType in ['scf', 'nscf', 'vc-relax']:
            occLocs = stri(evals, 'occupation numbers')

            #Check both numbers of k-point and occupation num. locations are 
            #equal
            assert len(kLocs) == len(occLocs), "The number of k-point and " + \
                    "'occupation number' occurences are not equal: {} and {}".\
                    format(len(kLocs), len(occLocs))
        else:
            occLocs = None

        return (kLocs, occLocs)


    def get_kpoint(self, evals, l):

        """
        Function returning the k-point vector in fractional
        coordinates given an index of its location in self.Lines

        Parameters:
        -----------
            kLoc: int
               Location of k-point information in self.Lines 

        Returns:
        --------
            kps: list of floats (3,)
                k-point components.
        """
        #Initialise 
        kps = []

        #String containing k-point information
        kp_string = evals[l].split("k =")[1].split("(")[0]
        #Extract the numbers and treat cases where numbers are joined
        #together
        for kval in kp_string.split():
            degeneracy = stri(kval, '-')
            if len(degeneracy) == 0:
                kps.append(float(kval))
            elif len(degeneracy) == 1:
                for dnum in degeneracy:
                    if dnum == 0:
                        kps.append(float(kval))
                    else:
                        kps.append(float(kval[:dnum]))
                        kps.append(float(kval[dnum:]))
            elif len(degeneracy) == 2:
                if degeneracy[0] == 0:
                    kps.append(float(kval[:degeneracy[1]]))
                    kps.append(float(kval[degeneracy[1]:]))
                else:
                    kps.append(float(kval[:degeneracy[0]]))
                    kps.append(float(kval[degeneracy[0]:degeneracy[1]]))
                    kps.append(float(kval[degeneracy[1]:]))
            elif len(degeneracy) == 3:
                kps.append(float(kval[:degeneracy[1]]))
                kps.append(float(kval[degeneracy[1]:degeneracy[2]]))
                kps.append(float(kval[degeneracy[2]:]))

        return kps


    def get_eigenvalues(self, kp_eigvals):

        """
        Returns the eigenvalues as a list when given as a list of strings where
        eigenvalues are space-separated, at a particular k-point.
        Parameters:
        -----------
            kp_eigvals: list of strings of space-separated eigenvalues
        Returns:
        --------
            eK: list of eigenvalues
        """

        eK = []
        for line in kp_eigvals:
            for val in line.split():
                degeneracy = stri(val, '-')
                if len(degeneracy) <= 1:
                    eK.append(float(val))
                elif len(degeneracy) == 2:
                    eK.append(float(val[:degeneracy[1]]))
                    eK.append(float(val[degeneracy[1]:]))
                else:
                    raise AssertionError("Too many values stuck together.")
        return eK


    def return_hyphen_separated_vals(self, val):
        """
        Function returns list of floats when they are separated by '-' in
        a string.
        """
        eK = []
        degeneracy = stri(val, '-')
        if degeneracy != []:
            for i, d in enumerate(degeneracy):
                if i + 1 != len(degeneracy):
                    if i == 0:
                        if d != 0:
                            eK.append(float(val[:d]))
                    eK.append(float(val[d:degeneracy[i+1]]))
                else:
                    if i == 0:
                        if d != 0:
                            eK.append(float(val[:d]))
                    eK.append(float(val[d:]))
        else:
            eK.append(float(val))
        return eK


    def get_qe_nonspinpol_bands(self):
        """
        Reads the bands in the non-spin polarised case
        """

        #Initialise dictionaries & temporary array
        e_bands_dict = {}
        k_points_dict = {}

        #Define section of 
        try:
            idx = strindex(self.Lines, ' k =', first=True)
            idx_end = strindex(self.Lines, 'the Fermi energy is')
            evals = self.Lines[idx + 2:idx_end - 1]
        except:
            AssertionError("Verbosity of nscf or bands calculation in" + \
                    " quantum-espresso must be set to 'high'")

        #Get k-point locations
        kLocs, occLocs = self.get_kpoint_locations(evals)

        for i, l in kLocs:
            #Get k-point
            k_points_dict[i+1] = self.get_kpoint(evals, l)

            #Get eigenvalues
            try:
                if self.calcType in ['scf', 'nscf', 'vc-relax']:
                    kp_eigvals = evals[kLocs[i] + 2:occLocs[i] - 1]
                elif self.calcType == "bands":
                    kp_eigvals = evals[kLocs[i] + 2:kLocs[i+1] - 1]
            except IndexError:
                if sp == 0:
                    end = strindex(evals, '\n') #Last occurrence
                elif sp == 1:
                    end = strindex(evals, 'Writing') #Spin down case
                    kp_eigvals = evals[kLocs[i]+2:end]

            e_bands_dict[i+1] = self.get_eigenvalues(kp_eigvals)

        #Re-order eigenvalues and k-points and shift by Fermi Energy
        e_bands_dicts[sp] = collections.OrderedDict(sorted(e_bands_dicts[sp].items()))
        e_bands = np.array(list(e_bands_dicts[sp].values()))

        self.eBands1 = e_bands - self.eFermi

        #Order k-point coordinates by k-point
        k_point_dict = collections.OrderedDict(sorted(k_points_dict.items()))
        k_points = np.array(list(k_point_dict.values()))
        self.kPoints = k_points

        return None


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

        return None


    def plot_bs(self):
        """
        This function plots the band structure using the bands parsed by the
        read_bands function.
        """

        s = 5
        #Figure title
        plotname = "{}.pdf".format(self.seed)
        #Plot band structure
        plt.figure(1)
        plt.ylabel("$E - E_{Fermi}$ / eV")
        plt.plot(np.array(self.eBands1), "-k", linewidth=0.8)
        plt.plot(np.array(self.eBands2), "--r", linewidth=0.3)
        plt.axhline(y = 0)
        plt.xlim((0,len(self.eBands1 - 2)))
        plt.ylim((-7,7))

        # Plot bandgap
        plt.axhline(y=self.LU, color="k", linestyle="dashed", linewidth=0.8, alpha=0.8)
        plt.axhline(y=self.HO, color="k", linestyle="dashed", linewidth=0.8, alpha=0.8)

        #Tweak and save
        plt.tight_layout()
        plt.savefig(plotname)

        return None


    def plot_bs_special(self):
        s = 5
        #Figure title
        plotname = "{}.pdf".format(self.seed)
        kpoints = self.kPoints
        kpoints_values = [list(item) for item in list(kpoints.values())]
        special_points = self.Path
        special_points_values = [list(item) for item in list(special_points.values())]
        special_points_keys = list(special_points.keys())

        # find indices where keys are the same 
        special_indices = [i for i, item in enumerate(kpoints_values) if item in special_points_values]

        #Plot band structure
        fig, ax = plt.subplots()

        plt.plot(np.array(self.eBands1), "-k", linewidth=0.85)
        plt.plot(np.array(self.eBands2), "--r", linewidth=0.3)

        # labels
        ax.set_ylabel(r"$E - E_{\mathrm{Fermi}}$ / eV", fontsize=12)
        ax.set_xlim((0,special_indices[-1]))
        ax.set_ylim((-7.5,3))

        # Plot special k points 
        ax.set_xticks(special_indices)
        ax.set_xticklabels(special_points_keys, rotation=45, fontsize=11)

        # Plot vertical lines at these points 
        for x in special_indices:
            ax.axvline(x, color="k", linewidth=0.8, alpha=0.8)

        # Plot bandgap
        ax.axhline(self.LU, color="k", linestyle="dashed", linewidth=0.8, \
                alpha=0.8)
        ax.axhline(self.HO, color="k", linestyle="dashed", linewidth=0.8, \
                alpha=0.8)
        plt.tight_layout()
        plt.savefig(plotname)

        return None


if __name__ == "__main__":

    #Get the parsed seedname
    fname, kwargs = get_cmdline_options()

    #Get the band gap
    bs = BandStructure(fname, **kwargs)
    bs.read_bands()
    bs.get_kpath_indices()
    bs.band_gap()
    print("Band gap: {:.3f} eV".format(bs.bandGap))

    #print(bs.eBands1)
    #print(bs.eBands2)
    bs.seed = bs.fname.replace(".bands.out","")
    bs.plot_bs()

