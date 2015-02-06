/*
 * datatypes.h
 *
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */

#ifndef DATATYPESs_H_
#define DATATYPESs_H_
#include<string>

struct triple_int {
    int arg1, arg2, arg3;
}; // 3d vect

struct pair_double {
    double arg1, arg2;
}; // 2d vect

struct Particle_f {
	// for fast LJ computations
    double x_f;
    double y_f;
    double z_f;
    double eps_with_adsorbate; // epsilon with adsorbate, defined by mixing rule
    double sig_with_adsorbate_squared; // sigma with adsorbate, defined by mixing rule
   // double reduced_mass; // for Feynman Hibbs
};

struct particle_g { 
// for guests particles
	// Cartesian coords
	double x;
	double y;
	double z;
	// fractional coords
	double x_f;
	double y_f;
	double z_f;
	// adsorbate identity
	int type;
};

struct GridParameters {
	// These are extracted from the simulation.input file
    string forcefieldname;
    bool verbose; // 0 or 1 for verbose printing
    double grid_resolution;
    string gridoutputformat;
    double r_cutoff_squared; // cutoff for LJ potential (A), squared
    bool feynmanhibbs;
    double energy_threshold; // for void fraction calc
    double T; // temperature (K), only needed if Feynman Hibbs is true

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

    // these must be computed separately and read from a unit cell replication file
    int replication_factor_a; // for replicating unit cells given by .cssr
    int replication_factor_b;
    int replication_factor_c;

    // TODO why not use framework.t_matrix? It will not pass correctly to cuda kernal. not sure why
    double t_matrix[3][3];
};

struct GCMCParameters {
	// These are extracted from the simulation.input file
    string forcefieldname;
    bool pocketblocking; // 1: enabled
    bool verbose; // 0 or 1 for verbose printing
    bool debugmode; // print details of each move for debugging
    double r_cutoff_squared; // cutoff for LJ potential (A), squared
    int feynmanhibbs;
    double T; // temperature (K)
    double delta; // spatial step in moves
    int numtrials, samplefrequency, numequilibriumtrials; // sample every X MC moves, burn period equilibrium_trials
    double p_move; double p_exchange; double p_identity_change; // probability of a move and exchange with bath(delete/insert). or identity change for dual component only

    // These are given as arguments to binary gcmc
    string frameworkname;
    string adsorbate[2]; // label in FF for GCMC
    int numadsorbates; // number of adsorbates
    double fugacity[2];

    // This must be read from force field object
    double epsilon_matrix[2][2]; double sigma_squared_matrix[2][2]; // pure sigma, epsilon for adsorbate

    // this is read from the masses.def file
    double adsorbateMW[2];

    // these are computed internally, as they are dependent on framework
    int N_framework_atoms;

    // these must be computed separately and read from a unit cell replication file
    int replication_factor_a; // for replicating unit cells given by .cssr
    int replication_factor_b;
    int replication_factor_c;

    // TODO why not use framework.t_matrix? It will not pass correctly to cuda kernal. not sure why
    double t_matrix[3][3];
    double inv_t_matrix[3][3];
};

struct Grid_info {
    // for GCMC ease
    int N_x, N_y, N_z;
    double dx_f, dy_f, dz_f;
    int numtotalpts;
    int numpockets[2]; // if pocket blocking...
};

struct GCMC_stats
{
	// count trials
	int N_insertion_trials;
	int N_deletion_trials;
	int N_move_trials;
	int N_ID_swap_trials;
	// count accepted monte carlo trials 
	int N_insertions;
	int N_deletions;
	int N_moves;
	int N_ID_swaps;
	// count samples
	int N_samples;
	// average energies
	double guest_guest_energy_avg;
	double framework_guest_energy_avg;
	// average number of guests
	double N_g_avg[2];
	double N_g2_avg[2]; // squared
};

#endif /* DATATYPES_H_ */
