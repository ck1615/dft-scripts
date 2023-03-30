#!/usr/bin/env python3
import numpy as np
from numpy.linalg import norm
from glob import glob
from misctools import strindex
from ase.units import Hartree, Bohr
from ase.io import read, write
from misctools import strindices as stri

#Plotting
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#Parsing command line options
import sys
import getopt


class BandStructure:
    """
    This class contains all the methods and attributes for a band structure
    plot
    """

    def __init__(self, fname, figsize=6, ratio=0.9, ymin=-1.5, ymax=5,
                 units='ase'):

        # File contents
        with open(fname, 'r+') as f:
            self.Lines = f.readlines()

        #for line in self.Lines:
        #    print(line)

        self.cellname = fname.replace("bands", "cell")

        # self.tolerance:
        self.tol = 1e-10

        # Plot characteristics
        self.ratio = ratio
        self.figsize = figsize
        self.spin_down_colour = ':ro'
        self.spin_up_colour = '--bo'
        self.markersize = self.figsize / 8
        self.linewidth = self.markersize / 3
        self.xlim = (0, 1)
        self.xlabel = 'Wavevectors'
        self.ylabel = r'$E - E_{\mathrm{Fermi}}$ / eV'
        self.y_majorticks = 1.0
        self.y_minorticks = 0.5
        self.y_major_tick_formatter = '{x:.0f}'
        self.ylim = (ymin, ymax)

        # Band gap characteristics
        self.HO = None
        self.LU = None
        self.band_gap = None
        self.nocc = None

        # Bands and path related variables
        self.kpath_idx = []   # Indices from 0.0 to 1.0 of k-points
        self.path = None
        self.path_ticks = None
        self.labels = []
        self.fermi_energy = None

        # K-points and cell
        self.k_points = []
        self.cell = np.array([])

        # Change rc params
        plt.rcParams['axes.labelsize'] = 2*self.figsize
        plt.rcParams['xtick.bottom'] = False
        plt.rcParams['font.size'] = 2*self.figsize

        # Get data
        self.get_calc_data()
        self.get_highsym_data()


    def get_calc_data(self):
        """
        This function extracts the:
            - number of k-points,
            - spin-polarisation,
            - number of occupied states,
            - total number of states,
            - Fermi energy
        """
        # Number of spin components
        num_spin = int(float(self.Lines[strindex(self.Lines,
                   'Number of spin components')].split()[-1]))
        if num_spin == 1:
            self.spin_polarised = False
        else:
            self.spin_polarised = True

        # Number of k-points
        self.nkpts = int(float(self.Lines[strindex(self.Lines,
                    'Number of k-points')].split()[-1]))

        # Number of occupied site
        self.nocc = int(float(self.Lines[strindex(self.Lines,
                    'Number of electrons')].split()[-1]))
        self.fermi_energy= float(self.Lines[strindex(self.Lines,
                                 'Fermi energies')].split()[-1]) * Hartree
        # Number of bands
        self.nBands = int(float(self.Lines[strindex(self.Lines,
                   'Number of eigenvalues')].split()[-1]))

        # Get unit cell
        self.get_unit_cell()

        # Get k-points
        self.get_kpoints()
        self.get_bands()

    def get_unit_cell(self):
        """
        This function gets the unit cell matrix from the .bands file.
        """
        # Get index of 'Unit cell vectors'
        idx_cell = strindex(self.Lines, 'Unit cell vectors')

        # Construct matrix
        self.cell = np.array([
            [float(val) for val in self.Lines[idx_cell + i + 1].split()]
            for i in range(3)])

        return None

    def get_kpoints(self):

        # Get the locations of the k-point lines
        k_indices = stri(self.Lines, 'K-point')

        # Initiate k-points array
        self.k_points = np.zeros((self.nkpts, 3))
        # Initialise array of k-points indices
        self.k_points_idcs = []

        for idx in k_indices:
            k_idx = int(self.Lines[idx].split()[1])-1
            self.k_points_idcs.append(k_idx)
            self.k_points[k_idx, :] = [float(val) for val in
                                       self.Lines[idx].split()[2:-1]]
        # Matrix multiply with unit cell to get cartesian k-points
        self.k_points_cart = np.matmul(self.k_points, self.cell)

        #print(self.k_points, self.k_points_cart)
        return None

    def get_bands(self):
        '''
        This function gets the upper and lower eigenvalues
        '''
        # Initialise arrays
        self.eBands1 = np.zeros((self.nkpts, self.nBands))
        if self.spin_polarised:
            self.eBands2 = np.zeros((self.nkpts, self.nBands))

            # Get locations of upper and lower spin-components
            sc1_idcs = stri(self.Lines, 'Spin component 1')
            #sc2_idcs = stri(self.Lines, 'Spin component 2')

            for i, idx in enumerate(sc1_idcs):
                band_idx = self.k_points_idcs[i]
                self.eBands1[band_idx,:] = self.Lines[idx+1:idx+1+self.nBands]
                self.eBands2[band_idx,:] = self.Lines\
                        [idx + self.nBands + 2:idx + 2*self.nBands + 2]

            # Scale to eV
            self.eBands1 *= Hartree
            self.eBands1 -= self.fermi_energy
            if self.spin_polarised:
                self.eBands2 *= Hartree
                self.eBands2 -= self.fermi_energy


    def get_highsym_data(self):
        """
        Gets all data relative to the high-symmetry points in the Brillouin
        zone required to perform the plot.
        """
        self.get_kpath_indices()
        self.get_highsym_kpoints()
        self.get_highsym_ticks()

        # Get band gap
        self.get_band_gap()
        self.get_klocs_band_gap()

        return None

    def get_kpath_indices(self):
        """
        This function takes a set of k-points and associates to it a list of
        non-negative real numbers corresponding to the distance from the
        initial k-point on the path.
        """
        path = []
        for i, kpt in enumerate(self.k_points_cart):
            if i == 0:
                path.append(0.0)
            else:
                path.append(path[i-1] + norm(kpt - \
                            self.k_points_cart[i-1]))

        # Normalise list between 0.0 and 1.0
        self.kpath_idx = [(idx - path[0]) / (path[-1] - path[0]) for idx \
                          in path]

    def get_highsym_kpoints(self):
        """
        Gets the high-symmetry points used to perform the band structure
        calculation
        """
        #Open cell input file
        with open(self.cellname, 'r+') as f:
            self.input_lines = f.readlines()

        #Get start and end indices
        start_idx = strindex(self.input_lines, "%BLOCK SPECTRAL_KPOINT_PATH")
        end_idx = strindex(self.input_lines, "%ENDBLOCK SPECTRAL_KPOINT_PATH")

        #Extract path in crystal coordinates and symbols if present
        self.path = np.array([[float(l) for l in line.split()[:3]] for
            line in self.input_lines[start_idx + 1: end_idx] if \
                    line.split() != []])
        #Get symbols for high-symmetry points
        self.get_highsym_symbols(start_idx, end_idx)

    def get_highsym_symbols(self, start_idx, end_idx):
        """
        Gets the symbols of the high-symmetry points if present
        """

        self.labels = []

        for i, line in enumerate(self.input_lines[start_idx + 1:end_idx]):
            if line.split() == []:
                continue
            last_char = line.split()[-1]
            if '!' in last_char:
                #Modify symbol to make it suitable for plots
                init_symbol = last_char.strip('!')
                symbol = self.get_plottable_symbol(init_symbol)
                self.labels.append(symbol)
            else:
                self.labels.append(tuple(self.path[i]))

    def get_plottable_symbol(self, symbol):
        """
        This function takes the raw capitalised symbol and returns the
        Greek symbol in LaTeX symbols if required.
        """

        bands_input_symbols = {'G': r'$\Gamma$', "G'": r"$\Gamma$'",
                               'G"': r'$\Gamma$'
                               }
        if symbol in bands_input_symbols:
            return bands_input_symbols[symbol]
        else:
            return symbol

    def get_highsym_ticks(self):
        """
        This function gets the locations of the high-symmetry points along
        the x-axis of the plot.
        """

        #Define high-symmetry point ticks:
        self.path_ticks = np.zeros(len(self.path))

        #Ensure first and last high-symmetry point correspond with start
        #and end of k-point list
        init_diff = norm(self.path[0] - self.k_points[0])
        final_diff = norm(self.path[-1] - self.k_points[-1])

        #print("initial %1.3f and final %1.3f" % (init_diff, final_diff))
        #print(self.path[0],self.Structure.k_points['crystal'][0])
        #print(self.path[-1],self.Structure.k_points['crystal'][-1])

        assert init_diff < self.tol and final_diff < self.tol,\
        "Initial and final are not what is expected:" #+\
        #"initial %1.3f and final %1.3f" % (init_diff, final_diff)

        #Set the values of the first and last ticks
        self.path_ticks[0] = 0.0
        self.path_ticks[-1] = 1.0

        #Initial k-point index
        kpt_idx = 1
        #Iterate over non-extremal high-symmetry points
        for ip, p in enumerate(self.path[1:-1]):
            #Iterate over k-points
            for ik, k in enumerate(self.k_points):
                #Only consider k-points above index
                if ik < kpt_idx:
                    continue
                if norm(k-p) < self.tol:
                    kpt_idx = ik + 1 #Start at next k-point after match
                    self.path_ticks[ip + 1] = self.kpath_idx[ik]
                    break

    def get_band_gap(self):
        """
        This function computes the band gap explicitly from the k-point
        eigenvalues if necessary.
        """
        #Get highest occupied and lowest unoccupied levels
        #for spin-polarised case
        if self.spin_polarised:
            self.HO = max([
                    self.eBands1[:,self.nocc - 1].max(),
                    self.eBands2[:,self.nocc - 1].max()
                    ])
            self.LU= min([
                    self.eBands1[:,self.nocc - 1].min(),
                    self.eBands2[:,self.nocc - 1].min()
                    ])
            #Get band gap
            self.band_gap = max(self.LU - self.HO, 0.0)
        else:
            self.HO = self.eBands1[:,self.nocc - 1].max()
            self.LU = self.eBands1[:,self.nocc].min()
            #Get band gap
            self.band_gap = max(self.LU - self.HO, 0.0)

    def get_klocs_band_gap(self, tol=1e-4):
        """
        Get the indices of the k-points at which the HO & LU occur
        """
        #Get indices of HO k-points
        #Spin-polarised case
        if self.spin_polarised:
            self.kpt_idx_HO = np.append(np.where(
                abs(self.eBands1[:, self.nocc - 1] - self.HO) \
                        < tol
                ), np.where(abs(
                    self.eBands2[:, self.nocc - 1]) < tol))
            self.kpt_idx_LU = np.append(np.where(
                abs(self.eBands1[:, self.nocc] - self.LU) < tol/10
                ), np.where(abs(
                    self.eBands2[:, self.nocc]) < tol))
        #Spin-unpolarised case
        else:
            self.kpt_idx_HO = np.where(
                abs(self.eBands1[:, self.nocc - 1] - self.HO) < \
                        tol
                )[0]
            self.kpt_idx_LU = np.where(
                abs(self.eBands1[:, self.nocc] - self.LU) < tol
                )[0]

    def plot_band_structure(self, save_pdf=True):
        """
        This function plots the band structure
    """

        # Start plot
        fig, ax = plt.subplots(figsize=(self.figsize*self.ratio, self.figsize))
        # Spin polarised case
        if self.spin_polarised:
            ax.plot(
                    self.kpath_idx,
                    self.eBands1,
                    self.spin_up_colour, label='Spin up',
                    linewidth=self.linewidth,
                    markersize=self.markersize
                    )
            ax.plot(
                    self.kpath_idx,
                    self.eBands2,
                    self.spin_down_colour,
                    label='Spin down',
                    linewidth=self.linewidth,
                    markersize=self.markersize,
                    alpha=0.4
                    )
        else:
            ax.plot(
                    self.kpath_idx,
                    self.eBands1,
                    self.spin_up_colour,
                    linewidth=self.linewidth,
                    markersize=self.markersize
                    )

        # Set energy (y) axis quantities
        ax.yaxis.set_major_locator(MultipleLocator(self.y_majorticks))
        ax.yaxis.set_major_formatter(self.y_major_tick_formatter)
        ax.yaxis.set_minor_locator(MultipleLocator(self.y_minorticks))
        ax.set_ylabel(self.ylabel)
        ax.set_ylim(self.ylim)

        # Set high-symmetry point quantities
        ax.set_xlim((0.0, 1.0))
        ax.set_xticks(self.path_ticks)
        ax.set_xticklabels(
                self.labels,
                rotation=0,
                fontsize=self.figsize * 2
                )

        # Plot vertical lines at each high-symmetry point
        for i, t in enumerate(self.path_ticks):
            ax.axvline(x=t, c='k', linewidth=self.linewidth)

        # Plot horizontal line at the origin (Fermi energy)
        ax.axhline(y=0.0, c='k', linestyle='--', linewidth=self.linewidth)

        # Locations of the LU points
        lu_indices = []

        # Plot additions for insulating band structure
        if (self.band_gap is not None) and (self.band_gap > 0.05):
            # Coloured in section of gap
            ax.axhspan(self.HO, self.LU, alpha=0.3, color='green')
            # Positions of HO & LU k-points
            for ho_idx in self.kpt_idx_HO:
                ax.plot(self.kpath_idx[ho_idx], self.HO, 'ko',
                        ms=self.markersize * 4)

            for lu_idx in self.kpt_idx_LU:
                lu_loc = self.kpath_idx[lu_idx]
                if (lu_loc != 1.0):
                    lu_indices.append(lu_loc)

                ax.plot(lu_loc, self.LU, 'ko', ms=self.markersize*4)

            # If empty sequence
            if len(lu_indices) == 1:
                lu_idx = lu_indices[0]
            else:
                if 0.0 in lu_indices:
                    lu_indices.remove(0.0)
                lu_idx = min(lu_indices)

            # Double arrow indicating band gap
            if self.band_gap > 0.3:
                plt.arrow(
                        lu_idx,
                        self.HO,
                        0.0,
                        self.band_gap,
                        length_includes_head=True,
                        shape='full',
                        color='r',
                        head_width=0.01,
                        head_length=0.1,
                        )
                plt.arrow(
                        lu_idx,
                        self.LU,
                        0.0,
                        -self.band_gap,
                        length_includes_head=True,
                        shape='full',
                        color='r',
                        head_width=0.01,
                        head_length=0.08,
                        )
            # Positions of band gap mention
            x_bg = lu_idx + 0.04
            y_bg = 0.5 * (self.HO + self.LU)
            plt.text(x_bg, y_bg, "{:1.2f} eV".format(self.band_gap),
                     va='center', fontsize=1.6*self.figsize)

        # Save figure
        if save_pdf:
            fig.tight_layout()
            fig.savefig(self.cellname.replace(".cell",".pdf"))
            #fig.savefig("%s/%s" %
            #            (self.thdir ,self.Structure.xmlname.replace('xml', 'pdf')))

        return None

def main():
    """
    This main function reads a xml output file generated by PWscf (quantum-
    espresso software package) and outputs the band-structure
    """

    # Get filename and any optional arguments
    fname, kwargs = command_line_options()

    # Read XML data file and band structure stuff
    BS = BandStructure(fname, figsize=6, ratio=0.8, ymin=-1.5, ymax=20)
    BS.plot_band_structure(save_pdf=True)

    return None


def command_line_options():
    """
    This function parses the command line options to get the filename.
    """

    # Get filename
    try:
        xmlname = sys.argv[1]
    except IndexError:
        raise IndexError("No filename has been provided.")
    # Other arguments
    argv = sys.argv[2:]

    # Iterate through options
    kwargs = {}
    opts, args = getopt.getopt(argv, "s:w:f:")
    for opt, arg in opts:
        if opt in ['-s', '--save-pdf']:
            if arg in ['true', 'True', 'TRUE']:
                arg = True
            elif arg in ['false', 'False', 'False']:
                arg = False
            kwargs['save_fig'] = arg

        elif opt in ['-w', '--energy-window']:
            kwargs['ymin'] = float(arg[1])
            kwargs['ymax'] = float(arg[3])
        elif opt in ['-f', '--figure-size']:
            kwargs['figsize'] = float(arg)

    return xmlname, kwargs

if __name__ == "__main__":
    main()
