/*
 * Forcefield.h
 *   loads information from force field file
 *   Stores Lennard-Jones potentials of form 4 eps * ( (sig/r)^12 - (sig/r)^6 )
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
#ifndef FORCEFIELD_H_
#define FORCEFIELD_H_
#include<vector>
#include<string>


class Forcefield {
public:
    Forcefield(std::string forcefield, bool verbose=false);

    int numinteractions; // number of interactions defined
    std::string name;
    std::vector<double> epsilon; // LJ parameter units: K
    std::vector<double> sigma; // LJ parameter units: A
    std::vector<std::string> identity; // labels of atoms that correspond to framework and adsorbate objects
};

#endif /* FORCEFIELD_H_ */
