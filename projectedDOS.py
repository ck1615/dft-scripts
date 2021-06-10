#!/usr/bin/env python3
"""
This script contains functions for plotting the projected Density of States
from (i) OptaDOS post-processing files, (ii) q-e projwfc.x output
files and (iii) VASP DOSCAR/PROCAR files.


Command-line execution instructions:
    -p --program: optados, q-e, vasp
    -t --type: dos, pdos
    -s --spin-polarised: true, false

Usage:

(i)  python3 projectedDOS.py seed -p program -t plot_type
(ii) projectedDOS.py seed -p program -t plot_type [if permission is granted]

If VASP is used, seed is going to be the seed of the plotname

Note for usage (ii) permissions must be granted to run it as an executable,
e.g. chmod +x projectedDOS.py

"""
from misctools import strindex, strindices
import numpy as np
from glob import glob
import matplotlib.pyplot as plt


class DoS():

    def __init__(self, seed, program="q-e", plotType="pdos",\
            spinPol=False):
        self.seed = seed
        self.program = program.lower()
        self.type = plotType.lower()
        self.spinPolarised = spinPol
        self.mode = None
        self.labels = None

        return None


    def get_file_contents_labels(self):
        """
        This function extracts the contents of the required files to plot the
        DoS/pDoS
        """

        if self.program == "optados":
            #Get filename 
            if self.type == "pdos":
                fname = "{}.pdos.dat".format(self.seed)

                #Extract file contents
                with open(fname, "r+") as f:
                    self.fContents = f.readlines()
                #Get pDOS labels
                self.get_optados_pdos_labels()

            elif self.type == "dos":
                fname = "{}.adaptive.dat".format(self.seed)
                #Extract file contents
                with open(fname, "r+") as f:
                    self.fContents = f.readlines()

        elif self.program == "q-e":
            if self.type == "pdos":
                #Get files and pDOS component labels
                files = glob("{}.pdos_*.dat".format(self.seed))
                self.labels = [' '.join(f.split('.')[-2].split('_')[-2:]) for f in \
                        files]

                #Get content
                self.fContents = []
                for fname in files:
                    print(fname)
                    with open(fname, "r+") as f:
                        self.fContents.append(f.readlines())

            elif self.type == "dos":
                #Filename of total DOS
                fname = "{}.pdos_tot".format(self.seed)
                with open(fname, "r+") as f:
                    self.fContents = f.readlines()

        elif self.program == "vasp":
            ValueError("VASP output files are not yet supported.")
            #Get DOSCAR content
            with open("vasprun.xml", "r+") as f:
                self.fContents = f.readlines()
        return None


    def get_optados_pdos_labels(self):
        """
        This function gets the labels of the different pDOS components

        Parameters:
        -----------

        Returns:
        --------

        """

        #Get the start of the enumeration of delimiters
        startString = '#|                    Partial Density of States -- Projectors                 |\n'
        startIdx = strindex(self.fContents, startString)

        #Get the labels delimiter
        delimiter='#+----------------------------------------------------------------------------+\n'
        #Get the line indices corresponding to all delimiters except last
        delimIdxs = strindices(self.fContents, delimiter, nmin=startIdx)
        self.lastDelim = delimIdxs[-1]

        #Get labels
        labels=[]

        for idx in delimIdxs[:-1]:
            #Get line with label info excluding first and last elements
            labelLine = self.fContents[idx+3].split()[1:-1] 

            if labelLine[-1] == 'Up':
                labels.append("{} {}".format(labelLine[0], labelLine[2]))
            else:
                pass
        self.labels = labels


    def get_fermi_energy(self):
        """
        This function extracts the OptaDOS calculated Fermi Energy when
        making an OptaDOS DOS or pDOS plot.
        """

        if self.program == "optados":
            with open('{}.odo'.format(self.seed), 'r+') as f:
                odoLines = f.readlines()

            pattern = "Fermi energy (Adaptive broadening)"
            FermiIdx = strindex(odoLines, pattern)
            FermiLine = odoLines[FermiIdx].split()
            self.FermiEnergy = float(FermiLine[-5])

        elif self.program == "vasp":
            self.FermiEnergy = float(self.fContents[5].split()[3])

        elif self.program == "q-e":
            print("Getting Fermi Energy")
            try:
                print("Opening")
                with open('{}.nscf.out'.format(self.seed), 'r+') as f:
                    mainLines = f.readlines()
            except FileNotFoundError:
                try:
                    with open('{}.scf.out'.format(self.seed), 'r+') as f:
                        mainLines = f.readlines()
                except FileNotFoundError:
                    FileNotFoundError("Could not find scf or nscf PWscf" \
                            "file.")

            if self.spinPolarised:
                pattern = "the spin up/dw Fermi energies are"
            else:
                pattern = "the Fermi energy is"

            FermiIdx = strindex(mainLines, pattern)
            FermiLine = mainLines[FermiIdx].split()
            self.FermiEnergy = float(FermiLine[-2])

        return None


    def get_data(self):
        """
        This function extracts the data required to plot the (projected or
        full) density of states.

        """

        #Get file contents
        self.get_file_contents_labels()

        #Get Fermi Energy
        self.get_fermi_energy()

        #OptaDOS
        if self.program == 'optados':
            if self.type == 'pdos':

                #Section of file containing columns of data
                dataLines = self.fContents[self.lastDelim + 1:]
                #Put data into an array form
                self.data = np.array([[float(l) for l in line.\
                        split()] for line in dataLines])
            elif self.type == "dos":
                ValueError("OptaDOS full DOS not yet supported.")

        #Quantum Espresso
        elif self.program == 'q-e':
            if self.type == 'pdos':
                #Get initial data set
                self.data = np.array([[float(l) - self.FermiEnergy
                    if k != 2 else -float(l) for
                    k, l in enumerate(line.split())] for line in
                  self.fContents[0][1:]])

                #Append remaining data set
                for i, content in enumerate(self.fContents[1:]):
                    self.data = np.append(self.data, np.array([[(-1)**(idx)*float(l) \
                            for idx, l in enumerate(line.split()[1:])] for \
                            line in content[1:]]), axis=1)

            elif self.type == 'dos':
                self.data = np.array([[float(l) for l in line.split()] for \
                        line in self.fContents])

        return None

    def plot_pDOS(self):
        """
        This function plots the pDOS.
        """
        #Get data
        #self.get_data()

        ##Calculate number of pDOS components
        #nPDOS = len(self.labels)

        #from matplotlib import cycler
        ##Configure matplotlib
        #plt.style.use('seaborn-ticks')
        #plt.rc('lines', linewidth=0.8)
        #plt.rc('axes', prop_cycle=cycler(color=['tab:blue', 'tab:orange',\
        #    'tab:green', 'tab:red', 'tab:purple', 'tab:cyan', 'tab:pink',
        #    'tab:olive', 'tab:gray']))

        #plt.figure(figsize=(16/2,9/2))
        #plt.xlabel(r"$E - E_{\mathrm{Fermi}}$ / eV")
        #plt.ylabel("Density of States")
        #plt.minorticks_on()
        #plt.plot(self.data[:,0] - self.FermiEnergy, self.data[:,1:])
        #plt.xlim((-20,20))
        ##plt.ylim((-35,35))
        #plt.legend(self.labels)

        #plt.axvline(x=0, color="k", linewidth=0.1)
        #plt.tight_layout()
        #plt.savefig('{}.PDOS.pdf'.format(self.seed))

        #Get data
        self.get_data()

        #Number of pdos components
        nPDOS = len(self.labels)

        #Configure matplotlib
        plt.style.use('seaborn-ticks')
        plt.rc('lines', linewidth=0.8)

        #Figure properties
        plt.figure(figsize=(8, 4.5))
        plt.xlabel(r"$E - E_{\mathrm{Fermi}}$ / eV")
        plt.ylabel("Density of States")
        plt.minorticks_on()

        cmap = plt.get_cmap('jet')
        colours = cmap(np.linspace(0, 1.0, nPDOS))
        plt.plot(self.data[:,0], self.data[:,1:])
        plt.legend(self.labels)

        plt.ylim((0,100))
        plt.axvline(x=0, color="k", linewidth=0.1)
        plt.tight_layout()
        plt.savefig('{}_pdos.pdf'.format(self.seed))






        return None

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
    opts, args = getopt.getopt(argv, "p:t:s:")
    for opt, arg in opts:
        if opt in ['-p', '--program']:
            program = arg
            assert program in ["optados", "q-e", "vasp"],\
                    "Chosen program is not one of: optados, q-e"+\
                    " or VASP"
            kwargs['program'] = program
        elif opt in ['-t', '--type']:
            plotType = arg
            kwargs['plotType'] = plotType
        elif opt in ['-s', '--spin-polarised']:
            spinPol = arg
            if spinPol.lower() == "true":
                spinPol = True
            elif spinPol.lower() == "false":
                spinPol = False
            else:
                Warning("Spin-polarised option was incorrectly specified. " + \
                        "Assuming non spin-polarised DoS or projected DOS.")
                spinPol = False
            kwargs['spinPol'] = spinPol
        else:
            Warning("The program or the type of plot has not been specified."+\
                    "Default values will be used.")

    return seed, kwargs

def main():
    """
    Main program
    """
    #from glob import glob
    #Get the filename, program and plot type 
    seed, kwargs = get_cmdline_options()

    #Initialise class
    calcDOS = DoS(seed, **kwargs)
    calcDOS.plot_pDOS()


if __name__ == "__main__":
    main()



