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

    # Create RawCIFs
    if not os.path.isdir('SymmetrisedCIFs'):
        os.mkdir('SymmetrisedCIFs')

    if not os.path.isdir('RawCIFs'):
        os.mkdir('RawCIFs')

    # Create all unsymmetrised CIF files
    for fname in glob("*.vc-relax.out"):
        print(fname)
        cifname = fname.replace('.vc-relax.out', '.cif')
        if not os.path.exists('RawCIFs/{}'.format(cifname)):
            write(cifname, read(fname))
            shutil.copy(cifname, 'RawCIFs')
        else:
            continue

        if not os.path.exists('SymmetrisedCIFs/{}'.format(cifname)):
            # Make symmetrised CIFs
            findsym_wrap(cifname, print_cif=True, axeso='abc', axesm='ab(c)',
                         lattol=tol, postol=tol)
            time.sleep(2)
            shutil.move(cifname, 'SymmetrisedCIFs')

    # Get the modes data
    os.chdir('RawCIFs')

    # Move the reference structure
    try:
        element = sys.argv[1]
    except IndexError:
        IndexError('No element was given')

    assert element in ['Cu', 'Mg'], "Element has to be Cu or Mg."

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

def extract_modes_data(HSfile):

    assert ".cif" in HSfile, "This is not a .cif file."

    LSfiles = [f for f in glob("*.cif") if 'norumpling' not in f]
    LSfiles = LSfiles[:2]

    mode_dict = {}
    Iso = isodistort(HSfile, silent=False)
    for fname in LSfiles:
        print(fname)
        Iso.load_lowsym_structure(fname)
        Iso.get_mode_labels()
        gm1p = Iso.modevalues['GM1+'][:2]
        x3p = Iso.modevalues['X3+']
        print(gm1p, x3p)
        mode_dict[fname] = {'GM1+': gm1p, 'X3+': x3p}
    Iso.close()
    np.save('mode_dict.npy', mode_dict, allow_pickle="TRUE")


if __name__ == "__main__":
        main()
