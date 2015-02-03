/*
 * datatypes.h
 *
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */

#ifndef DATATYPES_H_
#define DATATYPES_H_

struct Particle_f {
	// for fast LJ computations
    double x_f;
    double y_f;
    double z_f;
    double eps_with_adsorbate; // epsilon with adsorbate, defined by mixing rule
    double sig_with_adsorbate; // sigma with adsorbate, defined by mixing rule
    double reduced_mass; // for Feynman Hibbs
};

struct GridParameters {
	// These are extracted from the simulation.input file
    string forcefieldname;
    bool pocketblocking; // 1: enabled
    bool verbose; // 0 or 1 for verbose printing
    double grid_resolution;
    string gridoutputformat;
    double r_cutoff; // cutoff for LJ potential (A)
    int feynmanhibbs;
    double energy_threshold; // for void fraction calc
    double T; // temperature (K)

    // These are given as arguments to binary writegrid
    string frameworkname;
    string adsorbate; // label in FF for GCMC

    // This must be read from force field object
    double epsilon_guest; double sigma_guest; // pure sigma, epsilon for adsorbate

    // this is read from the masses.def file
    double adsorbateMW;

    // these are computed internally, as they are dependent on framework
    int N_framework_atoms;
    int N_z,N_y,N_x; // energy grid size
    double dx_f,dy_f,dz_f;

    // these must be computed separately and read from a unit cell replication file
    int replication_factor_a; // for replicating unit cells given by .cssr
    int replication_factor_b;
    int replication_factor_c;
};

#endif /* DATATYPES_H_ */
