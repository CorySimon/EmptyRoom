#!/bin/bash

###
#   Compute adsorption isotherm at pressures in pressures.txt
###
adsorbate="N2"

grab_isosteric_heat=true

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
    
    echo -e "\tWriting the grid"
 #     ./writegrid $structure $adsorbate
    ./writegrid $structure N_in_N2
 #     ./writegrid $structure Coulomb

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
    
    # define filename for compiled isotherm data
    forcefield=$(grep "Forcefield" simulation.input | awk '{print $2}')
    isothermfilename=output_files/${structure}_${adsorbate}_${forcefield}_isotherm_${T}K.csv
    echo "Crystal density: $rho kg/m3" > $isothermfilename
    if [ "$grab_isosteric_heat" = true ]; then
        echo "Pressure(Pa),Loading(mol/m3),Qst(kJ/mol)" >> $isothermfilename
    else
        echo "Pressure(Pa),Loading(mol/m3)" >> $isothermfilename
    fi
     
    cat pressures.txt | while read P
    do  
        # grab info from outputfiles
        volumetricloadingline=$(grep "moles/m3" output_files/${structure}_${adsorbate}_${P}Pa_${T}K_gcmc.out)
        volumetricloadingline=($volumetricloadingline)
        if [ "$grab_isosteric_heat" = true ]; then
            isostericheat=$(grep "Isosteric" output_files/${structure}_${adsorbate}_${P}Pa_${T}K_gcmc.out | awk '{print $6}')
            echo "$P,${volumetricloadingline[3]},${isostericheat}" >> $isothermfilename
        else
            echo "$P,${volumetricloadingline[3]}" >> $isothermfilename
        fi
    done
    echo "    Isotherms available in $isothermfilename"
done
