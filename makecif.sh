#!/bin/bash

[[ ! -d "SymmetrisedCIFs" ]] && mkdir SymmetrisedCIFs
[[ ! -d "RawCIFs" ]] && mkdir RawCIFs

for i in *vc-relax.out
do
  if [ ! -f "SymmetrisedCIFs/${i/vc-relax.out/cif}" ]
  then
    qe2cif.py $i -s
    mv ${i/vc-relax.out/cif} SymmetrisedCIFs/.
  fi
  if [ ! -f "RawCIFs/${i/vc-relax.out/cif}" ]
  then
    qe2cif.py $i
    mv ${i/vc-relax.out/cif} RawCIFs/.
  fi
done

