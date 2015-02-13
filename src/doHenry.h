#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib> // for "exit"
#include<sstream> // string stream
#include<vector>
#include <limits>
#include <complex>
#include<cmath>
#include<assert.h>

#ifndef doHenry_H 
#define doHenry_H

#define min_r .0000000000000001  // don't want to divide by zero...
// keep it at 64 blocks and 256 threads for random number generator
#define NUM_THREADS 256
#define NUM_BLOCKS 64

inline __device__ void DeviceFractionalToCartesian(double t_matrix[][3],
        double x_f, double y_f, double z_f,
        double & x, double & y, double & z)
{ // compute Cartesian coordinates from fractional
    x = t_matrix[0][0] * x_f + t_matrix[0][1] * y_f + t_matrix[0][2] * z_f;
    y = t_matrix[1][0] * x_f + t_matrix[1][1] * y_f + t_matrix[1][2] * z_f;
    z = t_matrix[2][0] * x_f + t_matrix[2][1] * y_f + t_matrix[2][2] * z_f;
}

__global__ void doHenry(
    curandStateMtgp32 * state,
    HenryParameters parameters,
    HenryStats * statistics,
    FrameworkParticle * framework_atoms)
{
    // each thread writes to this vector stored in the block.
    __shared__ double temp_energy[NUM_THREADS]; // entry is E_i
    __shared__ double temp_boltz[NUM_THREADS]; // entry is exp(-beta*E_i)
//  __shared__ double temp_occupiable[NUM_THREADS]; // 1 if occupiable, 0 if inoccupiable
    
    // initialize statistics
    if (threadIdx.x==0) {
        // each block gets its statistics structure.
        statistics[blockIdx.x].weighted_energy_sum = 0.0; // \sum E_i exp(-beta*E_i)
        statistics[blockIdx.x].canonical_sum = 0.0; // // \sum exp(-beta*E_i)
        statistics[blockIdx.x].void_fraction = 0.0;
        statistics[blockIdx.x].N_insertions = 0;
    }
  
    for (int cycle = 0; cycle < parameters.numinsertionsperthread; cycle++) {
        //
        // Throw dart in unit cell to get insertion location
        //
        double x_f = curand_uniform_double(&state[blockIdx.x]);
        double y_f = curand_uniform_double(&state[blockIdx.x]);
        double z_f = curand_uniform_double(&state[blockIdx.x]);

        __syncthreads();
        
        // Cartesian coords of insertion point
        double x, y, z; 
        DeviceFractionalToCartesian(parameters.t_matrix,
            x_f, y_f, z_f,
            x, y, z);
        
        //
        // Compute guest-framework energy at this point
        //
        double energy = 0.0;
        double x_framework, y_framework, z_framework;
        for (int i = -parameters.replication_factor_a; i <= parameters.replication_factor_a; i++) { // direction of x and a
            for (int j = -parameters.replication_factor_b; j <= parameters.replication_factor_b; j++) { // direction of y and b
                for (int k = -parameters.replication_factor_c; k <= parameters.replication_factor_c; k++) { // direction of z and c

                    for (int f = 0; f < parameters.N_framework_atoms; f ++) {

                        // fractional coordinates of framework atom under consideration
                        double x_f_framework = framework_atoms[f].x_f + i;
                        double y_f_framework = framework_atoms[f].y_f + j;
                        double z_f_framework = framework_atoms[f].z_f + k;

                        DeviceFractionalToCartesian(parameters.t_matrix,
                                x_f_framework, y_f_framework, z_f_framework,
                                x_framework, y_framework, z_framework);

                        double dx = x - x_framework; // distances between framework and sample point
                        double dy = y - y_framework;
                        double dz = z - z_framework;

                        double r2 = dx*dx + dy*dy + dz*dz; //r squared

                        r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r; // min radius to prevent blow-up

                        if (r2 < parameters.r_cutoff_squared) {  // energy contribution is nonzero only if within cutoff distance
                            double sig_over_r_squared = framework_atoms[f].sig_with_adsorbate_squared / r2;
                            double sig_over_r_sixth = sig_over_r_squared * sig_over_r_squared * sig_over_r_squared;
                            energy += 4.0 * framework_atoms[f].eps_with_adsorbate * sig_over_r_sixth * (sig_over_r_sixth - 1.0); // energy of framework atom id with guest_atom
                         } // end "if within cutoff..."
                     } // end loop over framework atoms
                } // end loop over c-direction unit cell replication
             } // end loop over b-direction unit cell replication
         } // end loop over a-direction unit cell replication

        //
        //  COLLECT BLOCK STATISTICS
        //
        temp_energy[threadIdx.x] = energy;
        temp_boltz[threadIdx.x] = exp(-energy / parameters.T);

        // see if inaccessible at this point for void fraction calc
//      if (energy < parameters.vf_threshold)
//          temp_occupiable[threadIdx.x] = 1;
//      else
//          temp_occupiable[threadIdx.x] = 0;   
        
        __syncthreads();  // wait for all threads to be done for getting block stats
    
        if (threadIdx.x==0) {
            // get sums along block
            for (int thrd=0; thrd < NUM_THREADS; thrd++) {
                statistics[blockIdx.x].canonical_sum += temp_boltz[thrd]; // sum of boltzman factors  (used to calc K_H)
                statistics[blockIdx.x].weighted_energy_sum += temp_energy[thrd] * temp_boltz[thrd]; // sum of boltzman weighted energy
//              statistics[blockIdx.x].void_fraction += 1.0 * temp_occupiable[thrd] / (NUM_THREADS * parameters.numinsertionsperthread); 
            }
            statistics[blockIdx.x].N_insertions += NUM_THREADS; // update number of insertions performed
        }
        __syncthreads(); // race condtion. might update temp vars before this is finished. keep sync here!
    } // end parameters.numinsertionsperthread cycles of insertions for the threads

//  if (threadIdx.x==0)
//  {
//      printf("Block %d avg ensemble avg Energy (K) = %lf \n ",blockIdx.x,statistics[blockIdx.x].weighted_energy_sum/statistics[blockIdx.x].canonical_sum);
//      printf("Block %d avg void fraction = %lf \n ",blockIdx.x,statistics[blockIdx.x].void_fraction);
//      printf("Block %d avg Henry (mol/kg-m3) = %lf \n ",blockIdx.x,statistics[blockIdx.x].canonical_sum/(NUM_THREADS*parameters.numinsertionsperthread))/8.314/parameters.T;
//  }
}
#endif
