#!/usr/bin/env python3
"""
This module contains functions for band path interpolation
"""
import numpy as np
from numpy.linalg import inv, norm

def interp_path(cell, path, spacing=0.02):

    #Get reciprocal rell
    recip_cell = 2*np.pi*inv(cell)
    interp=[]

    #Define interpolated points
    for ip, p in enumerate(path[:-1,:]):
        interp.append(p)
        #Create vector in path segment direction
        path_dir_full = path[ip + 1] - p
        #Normalise according to spacing
        nPoints = int(round(norm(np.matmul(recip_cell, path_dir_full) / spacing)))
        path_dir = path_dir_full / nPoints

        for num in range(1, nPoints):
            interp.append(p + num*path_dir)

    return interp

def write_interp(interp, code='qe'):

    with open("kpoints.txt", "w+") as f:
        count=1
        for point in interp:
            for i,p in enumerate(point):
                if i < 2:
                    f.write(str(p) + " ")
                else:
                    f.write(str(p) + " " + str(count) + "\n")
            count += 1

    return None

if __name__ == "__main__":
    path = np.array([[0.0, 0.0, 0.0],
    [0.0,0.0,  0.5],
    [0.0, 0.5,  0.5],
    [0.0,0.5,  0.0],
    [0.0,  0.0,  0.0],
    [-0.5,  0.0,  0.0],
    [-0.5,  0.0,  0.5],
    [0.0,  0.5,  0.0],
    [0.5,  0.5,  0.0],
    [0.0,  0.0,  0.0],
    [-0.5,  0.5,  0.5]])

    cell = np.array([[2.341850000000000,  -1.711300000000000,   0.000000000000000],
 [2.341850000000000,   1.711300000000000,   0.000000000000000],
 [-0.850027426760860,   0.000000000000000,   5.057869394691238]])

    spacing = 0.05

    interp = interp_path(cell, path, spacing=spacing)
    write_interp(interp)




