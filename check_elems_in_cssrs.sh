#!/bin/bash

# --- #
# Checks that atoms in .cssr files present in structures in structurelist.txt are in the force field file.
# --- #

# get elements in cssrs. simultaenously check .cssrs exist
cat structurelist.txt | while read f
do
	cat data/structures/$f.cssr >> temp.txt
	if [ ! -e data/structures/$f.cssr ]; then
		echo "missing ../sim_data/cssrs/$f.cssr"
		exit
		fi
done
cat temp.txt | awk '{if(NF==14) print $2;}' | sort -g | uniq >> temp00.txt
rm temp.txt

#		Check all elements are in force field
ff_line=$(grep "Forcefield" simulation.input)
ff=( ${ff_line} )
cat data/forcefields/${ff[1]}.def | awk '{if(NF==3) print $1;}' | sort -g | uniq >> temp01.txt

missing_elems=$(comm -23 temp00.txt temp01.txt)
n_missing=$(echo ${missing_elems} | awk '{print length}')
rm temp00.txt
rm temp01.txt
if [ ${n_missing} != 0 ]; then
	echo -e "Elements in cssrs, not in ../sim_data/forcefields/${ff[1]} \n----\n"
	echo ${missing_elems}
	echo -e "----\n"
fi
