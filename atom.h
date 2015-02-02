/*
 * atom.h
 *   Stores fractional coordinates, atomic mass, and identity of an atom
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
using namespace std;


#ifndef ATOM_H_
#define ATOM_H_

struct Atom {
	double mass; // atomic mass (amu)
	string identity; // atom identity
	double x_f, y_f, z_f; // fractional coordinates
};


#endif /* ATOM_H_ */
