#!/usr/bin/env python3
from glob import glob
from isodistortfile import isodistort
import numpy as np
import sys
import os, shutil
from findsymfile import findsym_wrap
from ase.io import read, write
import time


tol=1e-3
def main():
    """
    This function gets the mode amplitudes of the GM1+ and X3+ modes
    for all cif files, relative to the high-temperature phase with
    no rumpling above.
    """

    # Create RawCIFs and SymmetrisedCIFs directories if not present
    if not os.path.isdir('SymmetrisedCIFs'):
        os.mkdir('SymmetrisedCIFs')

    if not os.path.isdir('RawCIFs'):
        os.mkdir('RawCIFs')

    # Define filenames
    filenames = glob("*.vc-relax.out")

    # Create all CIF files
    print('Creating CIF files')
    for fname in filenames:
        cifname = fname.replace('.vc-relax.out', '.cif')
        # Raw CIF files (ASE output)
        if not os.path.exists('RawCIFs/{}'.format(cifname)):
            write(cifname, read(fname))
            shutil.copy(cifname, 'RawCIFs')
        else:
            continue
        # Symmetrised CIF files (FINDSYM output)
        if not os.path.exists('SymmetrisedCIFs/{}'.format(cifname)):
            # Make symmetrised CIFs
            findsym_wrap(cifname, print_cif=True, axeso='abc', axesm='ab(c)',
                         lattol=tol, postol=tol)
            time.sleep(2)
            shutil.move(cifname, 'SymmetrisedCIFs')

    # Get the modes data
    os.chdir('RawCIFs')

    # Move the reference structure
    element = get_element(filenames[0])

    # Define directory and high-symmetry file for mode decomposition
    if element == "Cu":
        HSfile = 'La2CuO4_HTT.9.4.norumpling.cif'

        norump_direc = "/Users/christopherkeegan/OneDrive - Imperial" +\
                       ' College London/Documents/PhD_Research/phd-project/' +\
                       'Structures/La2CuO4/PBESOL/q-e/{}'.format(HSfile)

        if not os.path.exists(HSfile):
            shutil.copy(norump_direc, '.')

    elif element == "Mg":
        HSfile = 'La2MgO4_HTT.4.0.norumpling.cif'

        norump_direc = "/Users/christopherkeegan/OneDrive - Imperial" +\
                       ' College London/Documents/PhD_Research/phd-project/' +\
                       'Structures/La2MgO4/PBESOL/Conventional/CIF/{}'.\
                       format(HSfile)

        if not os.path.exists(HSfile):
            shutil.copy(norump_direc, '.')

    #Extract modes
    extract_modes_data(HSfile)

def get_element(fname):
    '''
    This function gets the B-site ion given the a sample filename
    '''

    return fname.split('La2')[-1].split('O4')[0]

def extract_modes_data(HSfile):

    assert ".cif" in HSfile, "This is not a .cif file."

    # Get list of low-symmetry files
    LSfiles = [f for f in glob("*.cif") if 'norumpling' not in f]

    # If present, get the mode_dict and add to it
    if os.path.exists('mode_dict.npy'):
        mode_dict = np.load('mode_dict.npy', allow_pickle='TRUE').item()
        LSfiles = [i for i in LSfiles if i not in mode_dict]
    else:
        mode_dict = {}

    Iso = isodistort(HSfile, silent=True)
    for fname in LSfiles:
        if fname not in mode_dict:
            try:
                Iso.load_lowsym_structure(fname)
                try:
                    Iso.get_mode_labels()
                except:
                    continue
                gm1p = Iso.modevalues['GM1+'][:2]
                try:
                    x3p = Iso.modevalues['X3+']
                except KeyError:
                    x3p = [0.0]
                print(fname, gm1p, x3p)
                mode_dict[fname] = {'GM1+': gm1p, 'X3+': x3p}
            except:
                continue
    Iso.close()
    np.save('mode_dict.npy', mode_dict, allow_pickle="TRUE")


if __name__ == "__main__":
        main()
