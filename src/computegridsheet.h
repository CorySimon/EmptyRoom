/*
 * computegridsheet.h
 *  CUDA code to compute grid sheet at x_f = constant.
 *  Created on: Feb 3, 2015
 *      Author: corymsimon
 */
#include <cuda.h>
#include<assert.h>
using namespace std;
#include<string>
#include<cstdlib> // for "exit"
#include<vector>
#include <limits>
#include <complex>


#ifndef COMPUTEGRIDSHEET_H_
#define COMPUTEGRIDSHEET_H_

#define min_r .000000000001  // don't want to divide by zero...

inline __device__ void DeviceFractionalToCartesian(double t_matrix[][3],
		double x_f, double y_f, double z_f,
		double & x, double & y, double & z)
{ // compute Cartesian coordinates from fractional
    x = t_matrix[0][0] * x_f + t_matrix[0][1] * y_f + t_matrix[0][2] * z_f;
    y = t_matrix[1][0] * x_f + t_matrix[1][1] * y_f + t_matrix[1][2] * z_f;
    z = t_matrix[2][0] * x_f + t_matrix[2][1] * y_f + t_matrix[2][2] * z_f;
}

__global__ void ComputeGridSheet(
     double * z_f_gridpoints,
     double * y_f_gridpoints,
     double * zy_energies,
     FrameworkParticle * framework_atoms,
     GridParameters parameters,
     double x_f)
{
    double energy = 0.0; // each thread computes an energy at a particular point
    int z_index = threadIdx.x + blockIdx.x * blockDim.x; // which z point are we working on?
    int y_index = threadIdx.y + blockIdx.y * blockDim.y; // which y point are we working on?

    // if this thread does not have an assigned z or y point, do nothing.
    if ((z_index > (parameters.N_z - 1)) || (y_index > (parameters.N_y - 1)))
        return;

    // which fractional coordinate z_f and y_f is this thread responsible for?
    double z_f = z_f_gridpoints[z_index]; // local copy of z_f
    double y_f = y_f_gridpoints[y_index]; // local copy of y_f

    // set up Carteisan coordinates for insertion point and framework
    double x = 0.0, y = 0.0, z = 0.0; // Cartesian coords of grid point
    DeviceFractionalToCartesian(parameters.t_matrix, x_f, y_f, z_f, x, y, z);
    double x_framework = 0.0, y_framework = 0.0, z_framework = 0.0; // cartesian coords of framework (changes inside loop)

    for (int i = -parameters.replication_factor_a; i <=parameters.replication_factor_a; i++) { // direction of x and a
        for (int j = -parameters.replication_factor_b; j <= parameters.replication_factor_b; j++) { // direction of y and b
            for (int k = -parameters.replication_factor_c; k <= parameters.replication_factor_c; k++) { // direction of z and c

            	for (int framework_atom_index = 0; framework_atom_index < parameters.N_framework_atoms; framework_atom_index ++) {

                    // fractional coordinates of framework atom under consideration
                    double x_f_framework = framework_atoms[framework_atom_index].x_f + i;
                    double y_f_framework = framework_atoms[framework_atom_index].y_f + j;
                    double z_f_framework = framework_atoms[framework_atom_index].z_f + k;

                    DeviceFractionalToCartesian(parameters.t_matrix,
                    		x_f_framework, y_f_framework, z_f_framework,
                    		x_framework, y_framework, z_framework);

                    double dx = x - x_framework; // distances between framework and sample point
                    double dy = y - y_framework;
                    double dz = z - z_framework;

                    double r2 = dx*dx + dy*dy + dz*dz; //r squared

                    r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r; // min radius to prevent blow-up

                    if (r2 < parameters.r_cutoff_squared) {  // energy contribution is nonzero only if within cutoff distance
                        double sig_over_r_squared = framework_atoms[framework_atom_index].sig_with_adsorbate_squared / r2;
                        double sig_over_r_sixth = sig_over_r_squared * sig_over_r_squared * sig_over_r_squared;
                        energy += 4.0 * framework_atoms[framework_atom_index].eps_with_adsorbate * sig_over_r_sixth * (sig_over_r_sixth - 1.0); // energy of framework atom id with guest_atom
//                         if (parameters.feynmanhibbs == 1)
//                         {
//                                 energy += 4.0 * framework_atoms[framework_atom_index].epsilon * 48.508979 / 24.0 / framework_atoms[framework_atom_index].reduced_mass * sigma_over_r_sixth * (132.0 * sigma_over_r_sixth - 30.0 ) / r2 / parameters.T;
//                                 // 48.5 = hbar^2 /kb in A^2 units
//                         }
                     } // end "if within cutoff..."
                 } // end loop over framework atoms
            } // end loop over c-direction unit cell replication
         } // end loop over b-direction unit cell replication
     } // end loop over a-direction unit cell replication

     if ((z_index < parameters.N_z) && (y_index < parameters.N_y)) {  // write energies
         int energy_index_here = z_index + y_index * parameters.N_z; // WTF why do I need to do this instead of call directly?
         zy_energies[energy_index_here] = energy;
     }
}
#endif /* COMPUTEGRIDSHEET_H_ */
