#!/bin/bash

## Grabs Henry coefficients from simulations

adsorbate="Xe"

if [ -e Henry_results_${adsorbate}.csv ]; then
    rm Henry_results_${adsorbate}.csv
    fi

echo "Structure,KH(mol/kgPa),energy(K),crystal_density(kg/m3)" >> Henry_results_${adsorbate}.csv

cat structurelist.txt | while read f
do
    E_line=$(grep "<energy> (K)" output_files/${f}_${adsorbate}_henry.out)
    rho_line=$(grep "Crystal density" output_files/${f}_${adsorbate}_henry.out)
    kh_line=$(grep "Henry coefficient (mol/(kg-Pa)):" output_files/${f}_${adsorbate}_henry.out)
    
    # turn to vector-like
    kh_line=($kh_line)
    rho_line=($rho_line)
    E_line=($E_line)

    echo "$f,${kh_line[3]},${E_line[2]},${rho_line[3]}" >> Henry_results_${adsorbate}.csv
done
