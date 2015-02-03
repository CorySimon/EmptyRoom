/*
 * Forcefield.h
 *   loads information from force field file
 *   Stores Lennard-Jones potentials of form 4 eps * ( (sig/r)^12 - (sig/r)^6 )
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
using namespace std;
#include<vector>
#include<string>

#ifndef FORCEFIELD_H_
#define FORCEFIELD_H_

class Forcefield {
public:
	Forcefield(string forcefield, bool verbose=true);

    int nointeractions; // number of interactions defined
    string name;
    vector<double> epsilon; // LJ parameter units: K
    vector<double> sigma; // LJ parameter units: A
    vector<string> identity; // labels of atoms that correspond to framework and adsorbate objects
};

#endif /* FORCEFIELD_H_ */
