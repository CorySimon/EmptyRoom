/*
 * computegridsheet.h
 *  CUDA code to compute grid sheet at x_f = constant.
 *  -- van der waals grids
 *  -- electrostatic potential energy (Coulomb) grids
 *  Created on: Feb 3, 2015
 *      Author: corymsimon
 */
#include <cuda.h>
#include<assert.h>
#include<string>
#include<cstdlib>  // for "exit"
#include<vector>
#include <limits>
#include <complex>


#ifndef SRC_COMPUTEGRIDSHEET_H_
#define SRC_COMPUTEGRIDSHEET_H_

#define min_r .000000000001  // don't want to divide by zero...

inline __device__ void DeviceFractionalToCartesian(double t_matrix[][3],
        double x_f, double y_f, double z_f,
        double & x, double & y, double & z) {
    // compute Cartesian coordinates from fractional
    x = t_matrix[0][0] * x_f + t_matrix[0][1] * y_f + t_matrix[0][2] * z_f;
    y = t_matrix[1][0] * x_f + t_matrix[1][1] * y_f + t_matrix[1][2] * z_f;
    z = t_matrix[2][0] * x_f + t_matrix[2][1] * y_f + t_matrix[2][2] * z_f;
}

__global__ void ComputevdWGridSheet(
     double * z_f_gridpoints,
     double * y_f_gridpoints,
     double * zy_energies,
     FrameworkParticle * framework_atoms,
     GridParameters parameters,
     double x_f) {
    // compute energy on a sheet of grid points.
    // each thread computes an energy at a particular point
    double energy = 0.0;
    
    // which grid point is this thread working on?
    int z_index = threadIdx.x + blockIdx.x * blockDim.x;  // which z point are we working on?
    int y_index = threadIdx.y + blockIdx.y * blockDim.y;  // which y point are we working on?

    // if this thread does not have an assigned z or y point, do nothing.
    if ((z_index > (parameters.N_z - 1)) || (y_index > (parameters.N_y - 1)))
        return;

    // which fractional coordinate z_f and y_f is this thread responsible for?
    double z_f = z_f_gridpoints[z_index];  // local copy of z_f
    double y_f = y_f_gridpoints[y_index];  // local copy of y_f

    // set up Carteisan coordinates for this grid point
    double x = 0.0, y = 0.0, z = 0.0; // Cartesian coords of grid point
    DeviceFractionalToCartesian(parameters.t_matrix, x_f, y_f, z_f, x, y, z);
    
    for (int i = -parameters.replication_factor_a; i <=parameters.replication_factor_a; i++) { // direction of x and a
        for (int j = -parameters.replication_factor_b; j <= parameters.replication_factor_b; j++) { // direction of y and b
            for (int k = -parameters.replication_factor_c; k <= parameters.replication_factor_c; k++) { // direction of z and c

                for (int framework_atom_index = 0; framework_atom_index < parameters.N_framework_atoms; framework_atom_index ++) {

                    // fractional coordinates of framework atom under consideration
                    double x_f_framework = framework_atoms[framework_atom_index].x_f + 1.0 * i;
                    double y_f_framework = framework_atoms[framework_atom_index].y_f + 1.0 * j;
                    double z_f_framework = framework_atoms[framework_atom_index].z_f + 1.0 * k;
                    
                    // cartesian coords of framework 
                    double x_framework = 0.0, y_framework = 0.0, z_framework = 0.0; 
                    DeviceFractionalToCartesian(parameters.t_matrix,
                            x_f_framework, y_f_framework, z_f_framework,
                            x_framework, y_framework, z_framework);

                    double dx = x - x_framework; // distances between framework and grid point
                    double dy = y - y_framework;
                    double dz = z - z_framework;

                    double r2 = dx*dx + dy*dy + dz*dz; //r squared
                    r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r; // min radius to prevent blow-up

                    //
                    // Van der Waals energy
                    //
                    if (r2 < parameters.r_cutoff_squared) {  // energy contribution is nonzero only if within cutoff distance
                        double sig_over_r_squared = framework_atoms[framework_atom_index].sig_with_adsorbate_squared / r2;
                        double sig_over_r_sixth = sig_over_r_squared * sig_over_r_squared * sig_over_r_squared;
                        energy += 4.0 * framework_atoms[framework_atom_index].eps_with_adsorbate * sig_over_r_sixth * (sig_over_r_sixth - 1.0); // energy of framework atom id with guest_atom
//                         if (parameters.feynmanhibbs == 1) {
//                                 energy += 4.0 * framework_atoms[framework_atom_index].epsilon * 48.508979 / 24.0 / 
//                                           framework_atoms[framework_atom_index].reduced_mass * sigma_over_r_sixth * 
//                                           (132.0 * sigma_over_r_sixth - 30.0 ) / r2 / parameters.T;
//                                 // 48.5 = hbar^2 /kb in A^2 units
//                         }
                     } // end "if within cutoff..."
                 } // end loop over framework atoms
            } // end loop over c-direction unit cell replication
         } // end loop over b-direction unit cell replication
     } // end loop over a-direction unit cell replication
     
    //
    // Write energies to array
    //
    if ((z_index < parameters.N_z) && (y_index < parameters.N_y)) {
        int energy_index_here = z_index + y_index * parameters.N_z; // WTF why do I need to do this instead of call directly?
        //         printf("(Thread.x, Thread.y)=(%d,%d). (Block.x, Block.y)=(%d,%d). Responsible for y_f=%f, z_f=%f.E[%d]=%f\n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, y_f, z_f, energy_index_here, energy);
        zy_energies[energy_index_here] = energy;
    }
}

__global__ void ComputeCoulombGridSheet(
     double * z_f_gridpoints,
     double * y_f_gridpoints,
     double * zy_energies,
     FrameworkParticle * framework_atoms,
     GridParameters parameters,
     EWaldParameters ew_params,
     double * b1, double * b2, double * b3,
     double x_f) {
    // compute energy on a sheet of grid points.
    // each thread computes an energy at a particular point
    double energy = 0.0;
    
    // which grid point is this thread working on?
    int z_index = threadIdx.x + blockIdx.x * blockDim.x;  // which z point are we working on?
    int y_index = threadIdx.y + blockIdx.y * blockDim.y;  // which y point are we working on?

    // if this thread does not have an assigned z or y point, do nothing.
    if ((z_index > (parameters.N_z - 1)) || (y_index > (parameters.N_y - 1)))
        return;

    // which fractional coordinate z_f and y_f is this thread responsible for?
    double z_f = z_f_gridpoints[z_index];  // local copy of z_f
    double y_f = y_f_gridpoints[y_index];  // local copy of y_f

    // set up Carteisan coordinates for this grid point
    double x = 0.0, y = 0.0, z = 0.0; // Cartesian coords of grid point
    DeviceFractionalToCartesian(parameters.t_matrix, x_f, y_f, z_f, x, y, z);

    //
    // Short-range Coulomb energy
    //
    double energy_Coulomb_sr = 0.0;  // sr = short range
    for (int i = -parameters.replication_factor_a; i <=parameters.replication_factor_a; i++) { // direction of x and a
        for (int j = -parameters.replication_factor_b; j <= parameters.replication_factor_b; j++) { // direction of y and b
            for (int k = -parameters.replication_factor_c; k <= parameters.replication_factor_c; k++) { // direction of z and c

                for (int framework_atom_index = 0; framework_atom_index < parameters.N_framework_atoms; framework_atom_index ++) {

                    // fractional coordinates of framework atom under consideration
                    double x_f_framework = framework_atoms[framework_atom_index].x_f + 1.0 * i;
                    double y_f_framework = framework_atoms[framework_atom_index].y_f + 1.0 * j;
                    double z_f_framework = framework_atoms[framework_atom_index].z_f + 1.0 * k;
                    
                    // cartesian coords of framework 
                    double x_framework = 0.0, y_framework = 0.0, z_framework = 0.0; 
                    DeviceFractionalToCartesian(parameters.t_matrix,
                            x_f_framework, y_f_framework, z_f_framework,
                            x_framework, y_framework, z_framework);

                    double dx = x - x_framework; // distances between framework and grid point
                    double dy = y - y_framework;
                    double dz = z - z_framework;

                    double r2 = dx*dx + dy*dy + dz*dz; //r squared
                    r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r; // min radius to prevent blow-up

                    if (r2 < ew_params.cutoff_squared) {
                        double r = sqrt(r2);
                        energy_Coulomb_sr += framework_atoms[framework_atom_index].charge / r * erfc(r * sqrt(ew_params.alpha)) / (4 * M_PI * ew_params.eps0);
                    }
                } // end loop over framework atoms
            } // end loop over c-direction unit cell replication
        } // end loop over b-direction unit cell replication
    } // end loop over a-direction unit cell replication
     
     //
     // Long-range Coulomb energy
     //
     double energy_Coulomb_lr = 0.0;  // lr = long range
     // sum ovr k vectors
     for (int kx = -ew_params.kx; kx <= ew_params.kx; kx ++) {
         for (int ky = -ew_params.ky; ky <= ew_params.ky; ky ++) {
             for (int kz = -ew_params.kz; kz <= ew_params.kz; kz ++) {
                // continue if zero vector
                if ((kx == 0) & (ky == 0) & (kz == 0))
                    continue;

                // reciprocal lattice vector k = (k0, k1, k2) we are looking at
                // kx: amnt of lattice vector x
                double k0 = kx * b1[0] + ky * b2[0] + kz * b3[0]; 
                double k1 = kx * b1[1] + ky * b2[1] + kz * b3[1]; 
                double k2 = kx * b1[2] + ky * b2[2] + kz * b3[2]; 
                double mag_k_squared = k0*k0 + k1*k1 + k2*k2;

                //
                // Compute Structural factor S(k), which has real and imaginary part
                //
                double S_real_part = 0.0;
                double S_im_part = 0.0;
                for (int framework_atom_index = 0; framework_atom_index < parameters.N_framework_atoms; framework_atom_index ++) {
                    // cartesian coords of framework atom
                    double x_framework = 0.0, y_framework = 0.0, z_framework = 0.0; 
                    DeviceFractionalToCartesian(parameters.t_matrix,
                            framework_atoms[framework_atom_index].x_f, 
                            framework_atoms[framework_atom_index].y_f,
                            framework_atoms[framework_atom_index].z_f,
                            x_framework, y_framework, z_framework);

                    // S(k) structure factor (2 \pi taken care of in reciprocal lattice vectors)
                    S_real_part += framework_atoms[framework_atom_index].charge * 
                        cos(k0 * x_framework + k1 * y_framework + k2 * z_framework);
                    S_im_part += framework_atoms[framework_atom_index].charge * 
                        sin(k0 * x_framework + k1 * y_framework + k2 * z_framework);
                }
                // for the point charge
                double S_real_part_this = cos(k0 * x + k1 * y + k2 * z);
                double S_im_part_this = sin(k0 * x + k1 * y + k2 * z);

                double mag_S_squared = S_real_part * S_real_part_this + S_im_part * S_im_part_this;

                // add contribution to long-range Coulomb potential
                energy_Coulomb_lr += exp(- mag_k_squared/ 4.0 / ew_params.alpha) / mag_k_squared * mag_S_squared;
            } // end kz loop
        } // end ky loop
    } // end kx loop
    energy_Coulomb_lr = energy_Coulomb_lr / parameters.volume_unitcell / ew_params.eps0;

    energy = energy_Coulomb_sr + energy_Coulomb_lr;
    // subtract off self-interaction energy later
    
    //
    // Write energies to array
    //
    if ((z_index < parameters.N_z) && (y_index < parameters.N_y)) {
        int energy_index_here = z_index + y_index * parameters.N_z; // WTF why do I need to do this instead of call directly?
        //         printf("(Thread.x, Thread.y)=(%d,%d). (Block.x, Block.y)=(%d,%d). Responsible for y_f=%f, z_f=%f.E[%d]=%f\n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, y_f, z_f, energy_index_here, energy);
        zy_energies[energy_index_here] = energy;
    }
}

#endif // SRC_COMPUTEGRIDSHEET_H_
