/*
 * write_to_outputfile.h
 *
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */

#ifndef WRITE_TO_OUTPUTFILE_H_
#define WRITE_TO_OUTPUTFILE_H_
#include<assert.h>


void WriteSettingsToOutputfile(FILE * outputfile,
        GCMCParameters parameters,
        Framework framework,
        Forcefield forcefield,
        GridInfo grid_info,
        GridInfo Coulomb_grid_info,
        EWaldParameters ew_params,
        std::vector<int> uc_reps,
        boo::matrix<double> t_matrix,
        boo::matrix<double> inv_t_matrix,
        std::vector<std::string> adsorbate,
        std::vector<double> adsorbateMW,
        std::map<int, std::string> int_to_beadlabel,
        std::vector<Adsorbate> adsorbatetemplates,
        boo::matrix<double> epsilon_matrix,
        boo::matrix<double> sigma_squared_matrix)
{
    //
    // Write date
    //
    time_t t = time(0); // get time now
    struct tm * now = localtime(&t);
    fprintf(outputfile, "date: %d-%d-%d\n\n", (now->tm_year + 1900), (now->tm_mon + 1), now->tm_mday);

    fprintf(outputfile, "Grand-canonical Monte Carlo simulation.\n\n");
    fprintf(outputfile, "Temperature = %f K\n", parameters.T);

    // adsorbate information
    fprintf(outputfile, "Number of adsorbates: %d\n", parameters.numadsorbates);
    for (int i = 0; i < parameters.numadsorbates; i ++) {
        fprintf(outputfile, "    Adsorbate %s\n", adsorbate[i].c_str());
        fprintf(outputfile, "        molecular weight = %f\n", adsorbateMW[i]);
        fprintf(outputfile, "        number of LJ beads: %d\n", adsorbatetemplates[i].nbeads);
        for (int b = 0; b < adsorbatetemplates[i].nbeads; b++) {
            fprintf(outputfile, "          %s, bead type ID %d\n", int_to_beadlabel[adsorbatetemplates[i].beadtypes[b]].c_str(), adsorbatetemplates[i].beadtypes[b]);
            fprintf(outputfile, "             x = (%f, %f, %f)\n", adsorbatetemplates[i].bead_xyz(0, b), adsorbatetemplates[i].bead_xyz(1, b), adsorbatetemplates[i].bead_xyz(2, b));
        }
        fprintf(outputfile, "        number of point charges: %d\n", adsorbatetemplates[i].ncharges);
        for (int c = 0; c < adsorbatetemplates[i].ncharges; c++)
            fprintf(outputfile, "             charge = %f, x = (%f, %f, %f)\n", adsorbatetemplates[i].charges[c], 
                        adsorbatetemplates[i].charge_xyz(0, c), adsorbatetemplates[i].charge_xyz(1, c), adsorbatetemplates[i].charge_xyz(2, c));
    }


    //
    // Structure information
    //
    fprintf(outputfile, "\nStructure: %s\n", framework.name.c_str());
    fprintf(outputfile, "    Crystal density (kg/m3): %f\n", framework.density);
    fprintf(outputfile, "    a = %f A; b = %f A; c = %f A\n    alpha = %f, beta = %f, gamma = %f\n",
            framework.a, framework.b, framework.c,
            framework.alpha / M_PI * 180.0, framework.beta / M_PI * 180.0, framework.gamma / M_PI * 180.0);
    fprintf(outputfile, "    Atoms in framework: %d\n", framework.noatoms);
    fprintf(outputfile, "    Atom    # in framework\n");
    for (int a = 0 ; a < forcefield.numinteractions; a++) { // for every atom in force field, look for it in framework
        int count_atom_a = 0;
        for (int k = 0; k < framework.noatoms; k++) {
            if (framework.identity[k] == forcefield.identity[a])
                count_atom_a += 1;
        }
        if (count_atom_a != 0)
            fprintf(outputfile, "    %s    %d\n", forcefield.identity[a].c_str(), count_atom_a);
    }
    fprintf(outputfile, "\n    Tranformation matrix:\n");
    for (int i=0; i<3; i++) {
        fprintf(outputfile, "    ");
        for (int j=0; j<3; j++){
            fprintf(outputfile, " %f " , t_matrix(i, j));
        }
        fprintf(outputfile, "\n");
    }
    fprintf(outputfile, "\n    Inverse tranformation matrix:\n");
    for (int i=0; i<3; i++) {
        fprintf(outputfile, "    ");
        for (int j=0; j<3; j++){
            fprintf(outputfile, " %f " , inv_t_matrix(i, j));
        }
        fprintf(outputfile, "\n");
    }

    //
    // Forcefield information
    //
    fprintf(outputfile, "\nForcefield: %s\n", forcefield.name.c_str());
    fprintf(outputfile, "    Unit cell replication factors: %d %d %d\n", 
                    uc_reps[0], uc_reps[1], uc_reps[2]);
    fprintf(outputfile, "    Pure bead interactions:\n");
    for (int i = 0; i < int_to_beadlabel.size(); i ++)
        fprintf(outputfile, "      %s-%s: epsilon = %f K, sigma = %f A\n", 
                           int_to_beadlabel[i].c_str(), int_to_beadlabel[i].c_str(), 
                           epsilon_matrix(i, i), sqrt(sigma_squared_matrix(i, i)));
    if (int_to_beadlabel.size() > 1) {
        fprintf(outputfile, "    Mixed bead interactions:\n");
        for (int bx = 0; bx < int_to_beadlabel.size(); bx++) {
            for (int by = 0; by < bx; by++) {
                fprintf(outputfile, "      %s-%s: epsilon = %f K, sigma = %f A\n", 
                           int_to_beadlabel[bx].c_str(), int_to_beadlabel[by].c_str(), 
                           epsilon_matrix(bx, by), sqrt(sigma_squared_matrix(bx, by)));

            assert(sigma_squared_matrix(bx, by) == sigma_squared_matrix(by, bx));  // symmetry
            assert(epsilon_matrix(bx, by) == epsilon_matrix(by, bx));  // symmetry
            }
        }
    }

    //
    // GCMC simulation settings
    //
    fprintf(outputfile, "\nGCMC settings\n");
    fprintf(outputfile, "    Delta step for moves (A) = %f\n", parameters.delta);
    fprintf(outputfile, "    Total # Markov chain moves: %d\n", parameters.numtrials);
    fprintf(outputfile, "    Markov chain move sample frequency:  %d\n", parameters.samplefrequency);
    fprintf(outputfile, "    MC moves for equilibrium: %d\n", parameters.numequilibriumtrials);
    fprintf(outputfile, "    Prob(particle exchange): %f\n", parameters.p_exchange);
    fprintf(outputfile, "    Prob(particle translation): %f\n", parameters.p_move);
    fprintf(outputfile, "    Prob(particle ID swap): %f\n", parameters.p_identity_change);
    fprintf(outputfile, "    Prob(regrow): %f\n", parameters.p_regrow);
    
    if (parameters.charged_adsorbate_flag) {
        fprintf(outputfile, "\nEwald summation parameters.\n");
        fprintf(outputfile, "    alpha convergence parameter: %f\n", ew_params.alpha);
        fprintf(outputfile, "    short-range cutoff (A): %f\n", sqrt(ew_params.cutoff_squared));
        fprintf(outputfile, "    k-space replications = (%d, %d, %d)\n", ew_params.kx, ew_params.kx, ew_params.kz);
        fprintf(outputfile, "    reciprocal lattice vectors:\n");
        fprintf(outputfile, "       bx = [%f, %f, %f]\n", ew_params.b1[0], ew_params.b1[1], ew_params.b1[2]);   
        fprintf(outputfile, "       by = [%f, %f, %f]\n", ew_params.b2[0], ew_params.b2[1], ew_params.b2[2]);   
        fprintf(outputfile, "       bz = [%f, %f, %f]\n", ew_params.b3[0], ew_params.b3[1], ew_params.b3[2]);   
    }

    //
    // Grid info
    //
    fprintf(outputfile, "\nvdW Energy grid: %d by %d by %d points. Total grid points = %d\n",
            grid_info.N_x, grid_info.N_y, grid_info.N_z, grid_info.numtotalpts);
    double da = framework.a / (grid_info.N_x - 1); // actual grid spacing
    double db = framework.b / (grid_info.N_y - 1);
    double dc = framework.c / (grid_info.N_z - 1);
    fprintf(outputfile, "    Cartesian grid spacing ~: da = %f A, db = %f A, dc = %f A\n", da, db, dc);
    fprintf(outputfile, "    Pocking blocking: %d\n", parameters.pocketblocking);
    if (parameters.pocketblocking) {
        for (int i = 0; i < parameters.numadsorbates; i ++)
            fprintf(outputfile, "    %s grid, number of blocked pockets = %d\n", adsorbate[i].c_str(), grid_info.numpockets[i]);
    }
    if (parameters.charged_adsorbate_flag) {
        fprintf(outputfile, "\nCoulomb energy grid: "
                            "%d by %d by %d grid points. Total grid pts = %d\n", 
                        Coulomb_grid_info.N_x, Coulomb_grid_info.N_y, Coulomb_grid_info.N_x, Coulomb_grid_info.numtotalpts);
        fprintf(outputfile, "    Cartesian grid spacing ~(dx, dy, dz) = (%f, %f, %f)\n",
                        framework.a / (Coulomb_grid_info.N_x - 1),
                        framework.b / (Coulomb_grid_info.N_y - 1),
                        framework.c / (Coulomb_grid_info.N_z - 1));
    }
}

#endif /* WRITE_TO_OUTPUTFILE_H_ */
