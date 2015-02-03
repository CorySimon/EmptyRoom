#!/bin/bash
#
#  Simulates all structures in structurelist.txt
#

# check that cssrs are present and elements are present in force field
./check_elems_in_cssrs.sh

# get cutoff radius to feed to Julia script to compute UC replications
rc=($(grep "CutoffRadius" simulation.input))
rc=${rc[1]}


cat structurelist.txt | while read f
    do

    echo "Structure $f"
    
    # write grid
    echo -e "\tWriting grid"
    julia print_uc.jl $f once # get uc replications for one rc
    ./writegrid $f CH4
    rm data/uc_replications/$f.uc # remove uc rep file

 #     # run GCMC
 #     echo -e "\tRunning GCMC"
 #     julia print_uc.jl $f twice # get uc replications for two rc
 # 	./grand_canonical $f 20100.0 80600.0
 #     rm ../sim_data/uc_replications/$f.uc # remove uc rep file

    done
