#!/bin/bash
#
#  Simulates all structures in structurelist.txt
#

# check that cssrs are present and elements are present in force field
./check_elems_in_cssrs.sh

# get cutoff radius to feed to Julia script to compute UC replications
rc=($(grep "CutoffRadius" simulation.input))
rc=${rc[1]}

# get unit cell replications
julia print_uc.jl

cat structurelist.txt | while read f
    do

    echo "Structure $f"
    
 #     # write grid
 #     echo -e "\tWriting grid"
 #     ./writegrid $f Xe
 #     ./writegrid $f Kr
 
    # henry
    echo -e "\tComputing Henry"
    ./henry $f Xe

    # run GCMC
 #     echo -e "\tRunning GCMC"
 # 	./gcmc $f Xe 20100.0 Kr 80600.0
 #     rm data/uc_replications/$f.uc # remove uc rep file

    done
