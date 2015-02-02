/*
 * Framework.h
 *   Includes information about crystal structure, constructed by reading a cssr file
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
using namespace std;
#include<string>
#include<vector>
#include "atom.h"

#ifndef FRAMEWORK_H_
#define FRAMEWORK_H_

class Framework {
public:
	Framework(string cssrfilename, bool verbose=false); // constructor

	double a, b, c; // unit cell dimensions
	double alpha, beta, gamma; // unit cell angles
	int noatoms; // number of atoms
	vector<Atom> atoms; // atoms of the framework
	double density; // crystal density
	double volume_unitcell; // volume of unit cell (A^3)
};

#endif /* FRAMEWORK_H_ */
