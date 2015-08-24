/*
 * readsimulationmassesfile.h
 *   reads simulation.input file, unit cell replication factor file, and function to grab adsorbate MW
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */
#ifndef READSIMULATIONmassesfile_H_
#define READSIMULATIONmassesfile_H_
#include "Forcefield.h"
#include<cstdlib>
#include<cmath>
#include<cstdlib>
#include<string>
#include <fstream>
#include <sstream>


// assign adsorbate MW to parameters attribute
double GetAdsorbateMW(std::string adsorbate) {
    // open masses.def file with atomic masses
    std::string massesfilename = "data/masses.def";
    std::ifstream massesfile(massesfilename.c_str());
    if (massesfile.fail()) {
        printf("Masses file data/masses.def not present!\n");
        exit(EXIT_FAILURE);
    }

    double adsorbatemolecularweight = -1.0;
    std::string line;
    getline(massesfile, line); //waste line
    std::string label; double mass;
    while (getline(massesfile, line)) {
        std::istringstream linestream(line);
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

void ReadSimulationInputFile(GridParameters & parameters) {
    // Read simulation.input into parameters.
    std::string siminputfilename = "simulation.input";
    std::ifstream simfile(siminputfilename.c_str());
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
    parameters.EWald_precision = 1e-6;
    parameters.T = -1;
    
    // grab from simulation.input
    std::string word;
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
        if (word=="EWaldPrecision")
            simfile >> parameters.EWald_precision;
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

void ReadSimulationInputFile(HenryParameters & parameters) {
    // Read simulation.input into parameters. 
    std::string siminputfilename = "simulation.input";
    std::ifstream simfile(siminputfilename.c_str());
    if (simfile.fail())
    {
        printf("Simulation.input failed to open!\n");
        exit(EXIT_FAILURE);
    }

    // Initialize
    parameters.forcefieldname = "None";
    parameters.verbose = false;
    parameters.r_cutoff_squared = 12.5 * 12.5;
    parameters.T = -1;
    parameters.numinsertionsperA3 = -1;
    
    // grab from simulation.input
    std::string word;
    while (simfile >> word) {
        if (word=="Temperature")
            simfile >> parameters.T;
        if (word=="CutoffRadius") {
            simfile >> parameters.r_cutoff_squared;
            parameters.r_cutoff_squared = parameters.r_cutoff_squared * parameters.r_cutoff_squared;
        }
        if (word=="Forcefield")
            simfile >> parameters.forcefieldname;
        if (word=="NumberOfInsertionsPerA3")
            simfile >> parameters.numinsertionsperA3;
        if (word=="verbosemode")
            simfile >> parameters.verbose;
    }

    // check for missing
    if (parameters.T < 0.0) {
        printf("Missing Temperature in simulation.input");
        exit(EXIT_FAILURE);
    }
    if (parameters.numinsertionsperA3 < 0.0) {
        printf("Missing NumberOfInsertionsPerA3 in simulation.input");
        exit(EXIT_FAILURE);
    }
    if (parameters.forcefieldname == "None") {
        printf("Missing Forcefield in simulation.input");
        exit(EXIT_FAILURE);
    }
}

void ReadSimulationInputFile(GCMCParameters & parameters) {
    // Read simulation.input into parameters. bool gcmc is true if this is gcmc
    std::string siminputfilename = "simulation.input";
    std::ifstream simfile(siminputfilename.c_str());
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
    parameters.writeadsorbatepositions = false;
    parameters.EWald_precision = 1e-6;

    // grab from simulation.input
    std::string word;
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
        if (word=="DebugMode")
            simfile >> parameters.debugmode;
        if (word=="EWaldPrecision")
            simfile >> parameters.EWald_precision;
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
        if (word=="TrialProbability(Regrow)")
            simfile >> parameters.p_regrow;
        if (word=="move_delta")
            simfile >> parameters.delta;
        if (word=="MakeAssertions")
            simfile >> parameters.makeassertions;
        if (word=="WriteAdsorbatePositions")
            simfile >> parameters.writeadsorbatepositions;
        if (word=="WritePositionFrequency")
            simfile >> parameters.writepositionfrequency;
        if (word=="NumSnapshots")
            simfile >> parameters.num_snapshots;
    }

    // check for missing
    if (parameters.T < 0.0) {
        printf("Missing Temperature in simulation.input");
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
    if (parameters.numequilibriumtrials > parameters.numtrials) {
        printf("Eq trials > total trials!\n");
        exit(EXIT_FAILURE);
    }
    if ( (parameters.p_identity_change + parameters.p_move + parameters.p_exchange + parameters.p_regrow > 1.00000001) |
         (parameters.p_identity_change + parameters.p_move + parameters.p_exchange + parameters.p_regrow < 0.99999999) ) {
        printf("MC trial probabilities do not add to 1.0...\n");
        printf("\tp(ID change) = %f\n", parameters.p_identity_change);
        printf("\tp(move) = %f\n", parameters.p_move);
        printf("\tp(exchange) = %f\n", parameters.p_exchange);
        printf("\tp(regrow) = %f\n", parameters.p_regrow);
        exit(EXIT_FAILURE);
    }
}

std::vector<int> ReadUnitCellReplicationFile(std::string frameworkname, std::string once_or_twice) {
    // read unit cell replication factors from unit cell rep file
    // once_or_twice indicates cell needs to be twice the size of rc or once
    std::vector<int> uc_dims(3);
    uc_dims[0] = -1; uc_dims[1] = -1; uc_dims[2] = -1;

    std::string uc_filename = "data/uc_replications/" + frameworkname + "_" + once_or_twice + ".uc";
    std::ifstream ucfile(uc_filename.c_str());
    if (ucfile.fail()) {
        printf("unit cell replication factor file not found in data/uc_replications/$frameworkname_%s.uc\n", once_or_twice.c_str());
        exit(EXIT_FAILURE);
    }

    if ( !(ucfile >> uc_dims[0])) {
        printf("Problem reading UC file");
        exit(EXIT_FAILURE);
    }
    if ( !(ucfile >> uc_dims[1])) {
        printf("Problem reading UC file");
        exit(EXIT_FAILURE);
    }
    if ( !(ucfile >> uc_dims[2])) {
        printf("Problem reading UC file");
        exit(EXIT_FAILURE);
    }

    if ((uc_dims[0] == -1) | (uc_dims[1] == -1) | (uc_dims[2] == -1))
        printf("Error, problem reading unit cell replication file");
    
    return uc_dims;
}

std::vector<double> GrabGuestForceFieldParams(Forcefield forcefield, std::string adsorbate) {
    // searches through Forcefield object to get LJ params for adsorbate
    std::vector<double> eps_sig(2); // [epsilon, sigma] vector

    if (adsorbate == "Coulomb") {
        eps_sig[0] = 0.0;
        eps_sig[1] = 0.0;
        return eps_sig;
    }

    bool found_adsorbate = false;
    for (int i = 0; i < forcefield.numinteractions; i++) {
        if (forcefield.identity[i] == adsorbate) {
            found_adsorbate = true;
            eps_sig[0] = forcefield.epsilon[i];
            eps_sig[1] = forcefield.sigma[i];
            break;
        }
    }

    if (! found_adsorbate) {
        printf("Could not find adsorbate %s in forcefield file %s\n", adsorbate.c_str(), forcefield.name.c_str());
        exit(EXIT_FAILURE);
    }
    return eps_sig;
}

#endif /* READSIMULATIONmassesfile_H_ */
