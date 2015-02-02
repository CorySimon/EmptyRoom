/*
 * Framework.cpp
 *
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
#include<fstream>
#include<cstdlib>
#include<cmath>
#include <sstream>
#include "Framework.h"

Framework::Framework(string cssrfilename, bool verbose /*=false*/) {
	//
	// Read cssr file to import framework information
	//
	ifstream cssr(cssrfilename.c_str());
	if (cssr.fail()) {
		printf("CSSR file %s failed to import!", cssrfilename.c_str());
		exit(EXIT_FAILURE);
	}

	string line;
	// get cell dimensions on first line
	getline(cssr,line);
	istringstream input(line);
	input >> a >> b >> c;

	// get cell angles, convert to radians
	getline(cssr,line);
	input.str(line); input.clear();
	input >> alpha >> beta >> gamma;
	alpha = alpha * M_PI / 180.0;
	beta = beta * M_PI / 180.0;
	gamma = gamma * M_PI / 180.0;

	// get no of atoms
	getline(cssr,line);
	input.str(line); input.clear();
	input >> noatoms;

	// read atoms
	atoms.resize(noatoms); // give size of atom vector
	getline(cssr,line); // waste a line
	for (int atomnum = 0; atomnum < noatoms; atomnum++) {
		getline(cssr,line);
		input.str(line); input.clear();
		int junk; string junk2;
		input >> junk >> atoms[atomnum].identity >> atoms[atomnum].x_f >> atoms[atomnum].y_f >> atoms[atomnum].z_f;
		// reflect fractional coordiantes in [0,1]
		atoms[atomnum].x_f = fmod(atoms[atomnum].x_f, 1.0);
		atoms[atomnum].y_f = fmod(atoms[atomnum].y_f, 1.0);
		atoms[atomnum].z_f = fmod(atoms[atomnum].z_f, 1.0);
	}

	//
	// Read masses.def file to import atom masses
	//
	string massesfilename = "../sim_data/masses.def";
	ifstream massesfile(massesfilename.c_str());
    if (massesfile.fail()) {
        printf("File %s not present.\n", massesfilename.c_str());
        exit(EXIT_FAILURE);
    }

    getline(massesfile, line); //waste line
    // count lines
    int n_masses = 0;
    while (getline(massesfile, line))
    {
        n_masses ++;
    }

    massesfile.clear(); massesfile.seekg(0, massesfile.beg); // go back to beginning of file
    getline(massesfile, line); //waste line
    // to be extracted from masses.def
    vector<double> masses(n_masses);
    vector<string> identities(n_masses);
    for (int i = 0; i < n_masses; i++) {
		getline(massesfile, line); // get next line
		input.str(line); input.clear();
		input >> identities[i] >> masses[i];
	}

    //
    // Set masses in framework file object
    //
    for (int k = 0; k < noatoms; k++) {  // loop over framework atoms
        bool found = false; // ensure each atom is found
        for (int a = 0; a < n_masses; a++) {  // loop over masses.def
            if (atoms[k].identity == identities[a]) {  // if frameworkatom matches pseudoatom
                found = true;
                atoms[k].mass = masses[a];
                break;
            }
        }
        if (! found) {
            printf("Atom %s not present in %s.\n", atoms[k].identity.c_str(), massesfilename.c_str());
            exit(EXIT_FAILURE);
        }
    }

    //
    // Compute density
    //
    double v_unit_piped = sqrt(1.0 - pow(cos(alpha), 2.0) - pow(cos(beta), 2) - pow(cos(gamma), 2) + 2 * cos(alpha) * cos(beta) * cos(gamma)); // volume of unit parallelpiped
    volume_unitcell = v_unit_piped * a * b * c;
    double mass_unitcell = 0.0;
	for (int k = 0; k < noatoms; k++) {
		mass_unitcell += atoms[k].mass;
	}
	density = mass_unitcell/ volume_unitcell * 1660.53892; // (kg/m3) conversion factor for amu/A3 ---> kg/m3

    //
    // If verbose, print off info
    //
	if (verbose) {
		printf("a = %f, b = %f, c = %f\n", a, b, c);
		printf("alpha = %f, beta = %f, gamma = %f\n", alpha*180/M_PI, beta*180/M_PI, gamma*180/M_PI);
		printf("no of atoms: %d\n", noatoms);
		printf("Atom x_f y_f z_f mass\n---------------\n");
		for (int i = 0; i < noatoms; i++) {
			printf("%s %f %f %f %f\n", atoms[i].identity.c_str(), atoms[i].x_f, atoms[i].y_f, atoms[i].z_f, atoms[i].mass);
		}
		printf("\nCrystal density (kg/m3) = %f\n", density);
	}
}

