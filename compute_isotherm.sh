#!/bin/bash

###
#   Compute adsorption isotherm at pressures in pressures.txt
###
adsorbate="CH4"

# check that cssrs are present and elements are present in force field
./check_elems_in_cssrs.sh

# get cutoff radius to feed to Julia script to compute UC replications
rc=($(grep "CutoffRadius" simulation.input))
rc=${rc[1]}

# get unit cell replications
julia print_uc.jl $rc

# get number of processors
nprocs=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | tail -1)
echo "Number of processors used: $nprocs"

pressures=( $( cat pressures.txt ) ) # load pressure array
npressures=$(cat pressures.txt | wc -l)
T=$(grep "Temperature" simulation.input)
T=( $T )
T=${T[1]}

cat structurelist.txt | while read structure
do
    echo -e "\nComputing $adsorbate adsorbation isotherm in $structure at T=$T K"
    
 #     echo -e "\tWriting the grid"
 #     ./writegrid $structure $adsorbate

    for ((i=0; i<$npressures; i++))
    do  
        P=${pressures[i]}
        echo -e "\tComputing loading at P = $P Pa"
        # compute loading at each pressure with GCMC simulation
        ./gcmc ${structure} ${adsorbate} ${P} &
        
        # wait for jobs, since processors are full 
        [[ $((((i+1))%nprocs)) -eq 0 ]] && wait
    done

    wait

    # get crystal density 
    rho=$(grep "Crystal density" output_files/${structure}_${adsorbate}_${pressures[0]}Pa_${T}K_gcmc.out)
    rho=($rho)
    rho=${rho[3]} # kg/m3
    echo "rho=$rho"
    echo "Crystal density: $rho kg/m3" > output_files/${structure}_${adsorbate}_isotherm_${T}K.csv
    echo "Pressure(Pa),Loading(mol/m3)" >> output_files/${structure}_${adsorbate}_isotherm_${T}K.csv
     
    cat pressures.txt | while read P
    do  
        # grab info from outputfiles
        volumetricloadingline=$(grep "moles/m3" output_files/${structure}_${adsorbate}_${P}Pa_${T}K_gcmc.out)
        volumetricloadingline=($volumetricloadingline)
        echo "$P,${volumetricloadingline[3]}" >> output_files/${structure}_${adsorbate}_isotherm_${T}K.csv
    done
    echo "    Isotherms available in output_files/${structure}_${adsorbate}_isotherm_${T}K.csv"
done
