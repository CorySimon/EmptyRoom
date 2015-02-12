/*
 * write_to_outputfile.h
 *
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */

#ifndef WRITE_TO_OUTPUTFILE_H_
#define WRITE_TO_OUTPUTFILE_H_

void WriteSettingsToOutputfile(FILE * outputfile,
        GridParameters parameters,
        Framework framework,
        Forcefield forcefield,
        FrameworkParticle * framework_atoms) {
    //
    // Write date
    //
    time_t t = time(0); // get time now
    struct tm * now = localtime(&t);
    fprintf(outputfile, "date: %d-%d-%d\n\n", (now->tm_year + 1900), (now->tm_mon + 1), now->tm_mday);

    // adsorbate information
    fprintf(outputfile, "Adsorbate: %s\n", parameters.adsorbate.c_str());
    fprintf(outputfile, "    Molecular weight = %f\n\n", parameters.adsorbateMW);

    //
    // Structure information
    //
    fprintf(outputfile, "Structure: %s\n", framework.name.c_str());
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

    //
    // Forcefield information
    //
    fprintf(outputfile, "\nForcefield: %s\n", forcefield.name.c_str());
    fprintf(outputfile, "    %s epsilon (K) = %f, sigma (A) = %f\n", parameters.adsorbate.c_str(), parameters.epsilon_guest, parameters.sigma_guest);
    fprintf(outputfile, "    For LJ cutoff %f A, unit cell replication factors: %d %d %d\n", sqrt(parameters.r_cutoff_squared),
                    parameters.replication_factor_a, parameters.replication_factor_b, parameters.replication_factor_c);
    fprintf(outputfile, "    %s-XX interactions present:\n    XX    epsilon(K)    sigma(A)\n", parameters.adsorbate.c_str());
    vector<string> printed_these;
    int N_interactions = 0; // keep track of which interaction we printed
    for (int i = 0; i < framework.noatoms; i++) {
        string atom_id = framework.identity[i];
        if (i == 0) {
            N_interactions++;
            printed_these.push_back(atom_id);
            fprintf(outputfile, "    %s    %f    %f\n", atom_id.c_str(), framework_atoms[i].eps_with_adsorbate, sqrt(framework_atoms[i].sig_with_adsorbate_squared));
        }
        else {
            // check if already printed
            bool found = false;
            for (int k = 0; k < N_interactions; k++) {
                if (printed_these[k] == atom_id) {
                    found = true;
                    break;
                }
            }
            // print if not already
            if (! found) {
                N_interactions++;
                printed_these.push_back(atom_id);
                fprintf(outputfile, "    %s    %f    %f\n", atom_id.c_str(), framework_atoms[i].eps_with_adsorbate, sqrt(framework_atoms[i].sig_with_adsorbate_squared));
            }
        }
    }

    //
    // Grid info
    //
    fprintf(outputfile, "\nGrid: %d by %d by %d points. Total grid points = %d\n",
            parameters.N_x, parameters.N_y, parameters.N_z,
            parameters.N_x * parameters.N_y * parameters.N_z);
    double da = framework.a / (parameters.N_x - 1); // actual grid spacing
    double db = framework.b / (parameters.N_y - 1);
    double dc = framework.c / (parameters.N_z - 1);
    fprintf(outputfile, "    Actual grid spacing: da = %f A, db = %f A, dc = %f A\n", da, db, dc);
}


#endif /* WRITE_TO_OUTPUTFILE_H_ */
