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

Framework::Framework(string structurename, bool verbose /*=false*/) {
	//
	// Read cssr file to import framework information
	//
	ifstream cssr((structurename + ".cssr").c_str());
	if (cssr.fail()) {
		printf("CSSR file %s failed to import!", (structurename + ".cssr").c_str());
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
	x_f.resize(noatoms); // give size to vectors
    y_f.resize(noatoms);
    z_f.resize(noatoms);
    identity.resize(noatoms);
    mass.resize(noatoms);

	getline(cssr,line); // waste a line
	for (int i = 0; i < noatoms; i++) {
		getline(cssr,line);
		input.str(line); input.clear();
		int junk; string junk2;
		input >> junk >> identity[i] >> x_f[i] >> y_f[i] >> z_f[i];
		// reflect fractional coordiantes in [0,1]
		x_f[i] = fmod(x_f[i], 1.0);
		y_f[i] = fmod(y_f[i], 1.0);
		z_f[i] = fmod(z_f[i], 1.0);
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
            if (identity[k] == identities[a]) {  // if frameworkatom matches pseudoatom
                found = true;
                mass[k] = masses[a];
                break;
            }
        }
        if (! found) {
            printf("Atom %s not present in %s.\n", identity[k].c_str(), massesfilename.c_str());
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
		mass_unitcell += mass[k];
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
			printf("%s %f %f %f %f\n", identity[i].c_str(), x_f[i], y_f[i], z_f[i], mass[i]);
		}
		printf("\nCrystal density (kg/m3) = %f\n", density);
	}
}

