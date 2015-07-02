/*
 * Framework.h
 *   Includes information about crystal structure, constructed by reading a cssr file
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
#ifndef FRAMEWORK_H_
#define FRAMEWORK_H_
#include<string>
#include<vector>


class Framework {
public:
    // constructor
    Framework(std::string structurname, bool verbose=false); 

    std::string name;  // name of structure
    double a, b, c;  // unit cell dimensions
    double alpha, beta, gamma;  // unit cell angles
    int noatoms;  // number of atoms

    // fractional coords of atoms in crystal structure
    std::vector<double> x_f;
    std::vector<double> y_f;
    std::vector<double> z_f;

    // corresponding identities (elements)
    std::vector<std::string> identity; // identity

    // corresponding atomic masses
    std::vector<double> mass;  

    // corresponding charges
    std::vector<double> charge;  
    double net_charge;

    double density; // crystal density
    double volume_unitcell; // volume of unit cell (A^3)

    // for transferring between fractional and Cartesian coords
    double t_matrix[3][3]; // transformation matrix from fractional to Cartesian
    double inv_t_matrix[3][3]; // transformation matrix from Cartesian to fractional
};

#endif /* FRAMEWORK_H_ */
