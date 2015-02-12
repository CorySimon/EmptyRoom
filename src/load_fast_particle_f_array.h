/*
 * load_fast_particle_f_array.h
 *    Loads a malloc'ed pointer array of Particle_f datatype (see datatypes.h).
 *    This is recasting the information in Framework and Forcefield objects into a format
 *    that expedites energy calculations when looping over all framework atoms in the structure.
 *  Created on: Feb 3, 2015
 *      Author: corymsimon
 */

#ifndef LOAD_FAST_PARTICLE_F_ARRAY_H_
#define LOAD_FAST_PARTICLE_F_ARRAY_H_

void LoadFastFrameworkParticleArray(FrameworkParticle * framework_atoms, Framework framework, Forcefield forcefield, double epsilon_guest, double sigma_guest) {
    for (int i = 0; i < framework.noatoms; i++) {
        // keep framework_atoms on device as fractional coordinates [0,1]^3 (makes implementing periodic BC easier)
        framework_atoms[i].x_f = framework.x_f[i];
        framework_atoms[i].y_f = framework.y_f[i];
        framework_atoms[i].z_f = framework.z_f[i];

        bool found_frameworkatom = false;
        double epsilon; double sigma; // for storing epsilon and sigma for this framework atom
        for (int k = 0; k < forcefield.numinteractions;k++) {
            if (forcefield.identity[k] == framework.identity[i]) {
                epsilon = forcefield.epsilon[k];
                sigma = forcefield.sigma[k];
                found_frameworkatom = true;
                break;
            }
        }
        if ( ! found_frameworkatom) {
            printf("Could not find framework atom %s in forcefield.\n", framework.identity[i].c_str());
            exit(EXIT_FAILURE);
        }
        // implement mixing rules
        framework_atoms[i].sig_with_adsorbate_squared = pow((sigma_guest + sigma) / 2.0, 2);
        framework_atoms[i].eps_with_adsorbate = sqrt(epsilon_guest * epsilon);
        // get reduced mass with adsorbate for Feynman Hibbs
        // framework_atoms[f].reduced_mass = adsorbateMW * framework.atoms[f].mass / (adsorbateMW + framework.atoms[f].mass); // units: amu
    }
}

#endif /* LOAD_FAST_PARTICLE_F_ARRAY_H_ */
