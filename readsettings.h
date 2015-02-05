/*
 * readsimulationmassesfile.h
 *   reads simulation.input file, unit cell replication factor file, and function to grab adsorbate MW
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
using namespace std;
#include "Forcefield.h"
#include<cstdlib>
#include<cmath>
#include<cstdlib>
#include<string>
#include <fstream>
#include <sstream>
#ifndef READSIMULATIONmassesfile_H_
#define READSIMULATIONmassesfile_H_

// assign adsorbate MW to parameters attribute
double get_adsorbate_MW(string adsorbate) {
    // open masses.def file with atomic masses
    string massesfilename = "data/masses.def";
    ifstream massesfile(massesfilename.c_str());
    if (massesfile.fail()) {
        printf("Masses file data/masses.def not present!\n");
        exit(EXIT_FAILURE);
    }

    double adsorbatemolecularweight = -1.0;
    string line;
    getline(massesfile, line); //waste line
    string label; double mass;
    while (getline(massesfile, line))
    {
        istringstream linestream(line);
        linestream >> label >> mass;
        if (label == adsorbate)
        {
            adsorbatemolecularweight = mass;
            break;
        }
    }
    if (adsorbatemolecularweight < 0.0)
    {
        printf("No molecular weight for adsorbate %s in data/masses.def\n", adsorbate.c_str());
        exit(EXIT_FAILURE);
    }
    massesfile.close();
    return adsorbatemolecularweight;
}

void readsimulationinputfile(GridParameters & parameters) {
	// Read simulation.input into parameters. bool gcmc is true if this is gcmc
	string siminputfilename = "simulation.input";
	ifstream simfile(siminputfilename.c_str());
	if (simfile.fail())
	{
		printf("Simulation.input failed to open!\n");
		exit(EXIT_FAILURE);
	}

	// Initialize
    parameters.forcefieldname = "None";
    parameters.verbose = false;
    parameters.grid_resolution = -1.0;
    parameters.gridoutputformat = "txt";
    parameters.r_cutoff_squared = 12.5 * 12.5;
    parameters.feynmanhibbs = false;
    parameters.energy_threshold = 0.0; // for void fraction calc
    parameters.T = -1;
    
    // grab from simulation.input
    string word;
    while (simfile >> word) {
        if (word=="Temperature")
            simfile >> parameters.T;
        if (word=="VoidFractionEnergyThreshold")
            simfile >> parameters.energy_threshold;
        if (word=="GridResolution")
            simfile >> parameters.grid_resolution;
        if (word=="CutoffRadius") {
            simfile >> parameters.r_cutoff_squared;
            parameters.r_cutoff_squared = parameters.r_cutoff_squared * parameters.r_cutoff_squared;
        }
        if (word=="Forcefield")
            simfile >> parameters.forcefieldname;
        if (word=="GridOutputFormat")
            simfile >> parameters.gridoutputformat;
        if (word=="FeynmanHibbs")
            simfile >> parameters.feynmanhibbs;
        if (word=="verbosemode")
            simfile >> parameters.verbose;
    }

    // check for missing
    if (parameters.feynmanhibbs & (parameters.T < 0.0)) {
        printf("Missing Temperature in simulation.input");
        exit(EXIT_FAILURE);
    }
    if (parameters.grid_resolution < 0.0) {
        printf("Missing GridResolution in simulation.input");
        exit(EXIT_FAILURE);
    }
    if (parameters.forcefieldname == "None") {
        printf("Missing Forcefield in simulation.input");
        exit(EXIT_FAILURE);
    }
}

void readsimulationinputfile(GCMCParameters & parameters) {
	// Read simulation.input into parameters. bool gcmc is true if this is gcmc
	string siminputfilename = "simulation.input";
	ifstream simfile(siminputfilename.c_str());
	if (simfile.fail())
	{
		printf("Simulation.input failed to open!\n");
		exit(EXIT_FAILURE);
	}

	// Initialize
    parameters.forcefieldname = "None";
    parameters.pocketblocking = false;
    parameters.verbose = false;
    parameters.r_cutoff_squared = 12.5 * 12.5;
    parameters.feynmanhibbs = 0;
    parameters.T = -1;
    parameters.delta = 0.1;
    parameters.numtrials = -1;
    parameters.samplefrequency = -1;
    parameters.numequilibriumtrials = -1;

    // grab from simulation.input
    string word;
    while (simfile >> word) {
        if (word=="Temperature")
            simfile >> parameters.T;
        if (word=="CutoffRadius") {
            simfile >> parameters.r_cutoff_squared;
            parameters.r_cutoff_squared = parameters.r_cutoff_squared * parameters.r_cutoff_squared;
        }
        if (word=="Forcefield")
            simfile >> parameters.forcefieldname;
        if (word=="FeynmanHibbs")
            simfile >> parameters.feynmanhibbs;
        if (word=="verbosemode")
            simfile >> parameters.verbose;
        if (word=="PocketBlocking")
            simfile >> parameters.pocketblocking;
        if (word=="NumberOfTrials")
            simfile >> parameters.numtrials;
        if (word=="SampleFrequency")
            simfile >> parameters.samplefrequency;
        if (word=="CyclesForEquilibrium")
            simfile >> parameters.numequilibriumtrials;
        if (word=="TrialProbability(Move)")
            simfile >> parameters.p_move;
        if (word=="TrialProbability(Exchange)")
            simfile >> parameters.p_exchange;
        if (word=="TrialProbability(IdentityChange)")
            simfile >> parameters.p_identity_change;
        if (word=="move_delta")
            simfile >> parameters.delta;
    }

    // check for missing
    if (parameters.T < 0.0) {
        printf("Missing Temperature in simulation.input");
        exit(EXIT_FAILURE);
    }
    if (parameters.grid_resolution < 0.0) {
        printf("Missing GridResolution in simulation.input");
        exit(EXIT_FAILURE);
    }
    if (parameters.forcefieldname == "None") {
        printf("Missing Forcefield in simulation.input");
        exit(EXIT_FAILURE);
    }
    if ((parameters.numtrials == -1) | (parameters.samplefrequency== -1) | (parameters.numequilibriumtrials == -1)) {
        printf("Missing samplefrequency, CyclesForEquilibrium, or NumberOfTrials in simulation.input\n");
        exit(EXIT_FAILURE);
    }
}

triple_int readunitcellreplicationfile(string frameworkname) {
	// read unit cell replication factors from unit cell rep file
    triple_int uc_dims;

	uc_dims.arg1 = -1; // initialize
    uc_dims.arg2 = -1;
	uc_dims.arg3 = -1;
	string uc_filename = "data/uc_replications/" + frameworkname + ".uc";
	ifstream ucfile(uc_filename.c_str());
	if (ucfile.fail())
	{
		printf("unit cell replication factor file not found in data/uc_replications/$frameworkname.uc\n");
		exit(EXIT_FAILURE);
	}

	if ( !(ucfile >> uc_dims.arg1))
	{
		printf("Problem reading UC file");
		exit(EXIT_FAILURE);
	}
	if ( !(ucfile >> uc_dims.arg2))
	{
		printf("Problem reading UC file");
		exit(EXIT_FAILURE);
	}
	if ( !(ucfile >> uc_dims.arg3))
	{
		printf("Problem reading UC file");
		exit(EXIT_FAILURE);
	}

	if (uc_dims.arg1 == -1) printf("Error, problem reading unit cell replication file");
    
    return uc_dims;
}

void get_guest_FF_params_from_Forcefield(Forcefield forcefield, string adsorbate, GridParameters & parameters) {
	// searches through Forcefield object to get LJ params for adsorbate
	bool found_adsorbate = false;
	for (int i =0; i < forcefield.nointeractions; i++) {
		if (forcefield.identity[i] == adsorbate) {
			found_adsorbate = true;
			parameters.epsilon_guest = forcefield.epsilon[i];
			parameters.sigma_guest = forcefield.sigma[i];
			break;
		}
	}

	if (! found_adsorbate) {
		printf("Could not find adsorbate %s in forcefield file %s", adsorbate.c_str(), forcefield.name.c_str());
		exit(EXIT_FAILURE);
	}
}
#endif /* READSIMULATIONmassesfile_H_ */
