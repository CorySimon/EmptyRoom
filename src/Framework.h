/*
 * Framework.h
 *   Includes information about crystal structure, constructed by reading a cssr file
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
#ifndef FRAMEWORK_H_
#define FRAMEWORK_H_
using namespace std;
#include<string>
#include<vector>


class Framework {
public:
	Framework(string structurname, bool verbose=false); // constructor

	string name; // name of structure
	double a, b, c; // unit cell dimensions
	double alpha, beta, gamma; // unit cell angles
	int noatoms; // number of atoms
    // store atoms in crystal structure
	vector<double> x_f; // fractional coord
	vector<double> y_f; // fractional coord
	vector<double> z_f; // fractional coord
    vector<double> mass;
    vector<string> identity; // identity

	double density; // crystal density
	double volume_unitcell; // volume of unit cell (A^3)

	// for transferring between fractional and Cartesian coords
    double t_matrix[3][3]; // transformation matrix from fractional to Cartesian
    double inv_t_matrix[3][3]; // transformation matrix from Cartesian to fractional
};

#endif /* FRAMEWORK_H_ */
