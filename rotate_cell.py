#!/usr/bin/env python3
"""
This function rotates by 45° a unit cell for a Quantum ESPRESSO calculation.
"""
from ase import Atoms
from glob import glob
import os
from ase.io import read, write


def parse_seed():
    """
    This function parses arguments input in the terminal straight after the
    module.

    Returns:
    --------
    seed: str
    The seedname of all files
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("seed", type=str)
    seed = parser.parse_args().seed
        
    return seed

def rotate_45(seed):
    """
    This function rotates by 45 degrees
    """

    calc_type = seed.split(".")[-2]
    Atoms = read(seed)
    Atoms.rotate([1,0,0], [1,1,0], rotate_cell=False)
    write(seed.replace(".{}.in".format(calc_type), "_rot45.{}.in".
        format(calc_type))


if __name__ == "__main__":
    seed = parse_seed()
    rotate_45(seed)
