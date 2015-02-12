/*
 * Forcefield.cpp
 *
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */

#include "Forcefield.h"
#include<vector>
#include<string>
#include<fstream>
#include<cstdlib>
#include <sstream>

Forcefield::Forcefield(string forcefield, bool verbose /*=false*/) {
	// load forcefield
    name = forcefield;
    ifstream forcefieldfile(("data/forcefields/" + forcefield + ".def").c_str());
    if (forcefieldfile.fail()) {
    	printf("Forcefield file %s did not load.\n", ("data/forcefields/" + forcefield + ".def").c_str());
        exit(EXIT_FAILURE);
    }

    // count interactions
    nointeractions = 0;
    string line;
    getline(forcefieldfile, line); // waste a line
    while (getline(forcefieldfile, line)) nointeractions++;

    forcefieldfile.clear(); forcefieldfile.seekg(0, forcefieldfile.beg); // go back to beginning of file
    getline(forcefieldfile,line); // waste a line

    // define vector sizes
    epsilon.resize(nointeractions);
    sigma.resize(nointeractions);
    identity.resize(nointeractions);

    // extract LJ parameters
    for (int k = 0; k < nointeractions; k++) {
        getline(forcefieldfile, line);
        istringstream stringstream(line);
        stringstream >> identity[k] >> epsilon[k] >> sigma[k]; /// CHECK ORDER OF SIG EPS
    }

    // if verbose, print result
    if (verbose) {
        printf("Atom epsilon(K) sigma(A)\n");
    	for (int k = 0; k < nointeractions; k++) {
        	printf("%s %f %f\n", identity[k].c_str(), epsilon[k], sigma[k]);
        }
    }

}

