using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include<cuda.h>
#include <sys/time.h>
#include <cuda_runtime.h>
#include "datatypes.h"
#include "readsettings.h"
#include "Framework.h"
#include "Forcefield.h"
#include "write_settings_to_outputfile.h"
#include "load_fast_particle_f_array.h"
#include <curand.h>
#include <curand_kernel.h>
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include "doHenry.h"


#define CUDA_CALL(x) do { cudaError_t error = x;  \
  if (error != cudaSuccess) {  \
  printf("Error at %s:%d - %s \n",__FILE__,__LINE__, cudaGetErrorString(error)); \
  return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
		    return EXIT_FAILURE;}} while(0)

double read_timer() {
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

int main(int argc, char *argv[])
{
    if (argc != 3) {
    	printf("Run as:\n./henry structure_name AdsorbateID\n");
    	exit(EXIT_FAILURE);
    }
	//
	//  Import settings
	//
    HenryParameters parameters;
    parameters.frameworkname = argv[1];
    parameters.adsorbate = argv[2];
    parameters.adsorbateMW = get_adsorbate_MW(parameters.adsorbate);

    readsimulationinputfile(parameters);
    if (parameters.verbose) printf("Read simulation.input\n");
    triple_int uc_reps = readunitcellreplicationfile(parameters.frameworkname, "once");
    parameters.replication_factor_a = uc_reps.arg1;
    parameters.replication_factor_b = uc_reps.arg2;
    parameters.replication_factor_c = uc_reps.arg3;
    if (parameters.verbose) printf("Read .uc replication file\n");

    //
    // Construct forcefield and framework objects
    //
    Forcefield forcefield(parameters.forcefieldname);
    if (parameters.verbose) printf("Constructed Forcefield object\n");
    Framework framework(parameters.frameworkname);
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++){
			parameters.t_matrix[i][j] = framework.t_matrix[i][j];
		}
	}// TODO: remove and use framework.t_matrix instead. right now it cant pass to cuda kernel without memory error...
    parameters.N_framework_atoms = framework.noatoms;
    if (parameters.verbose) printf("Constructed Framework object\n");

    parameters.numinsertions = (int) parameters.numinsertionsperA3 * framework.volume_unitcell;

    // grab sigma/epsilon of adsorbate
    pair_double eps_sig = get_guest_FF_params_from_Forcefield(forcefield, parameters.adsorbate);
    parameters.epsilon_guest = eps_sig.arg1;
    parameters.sigma_guest = eps_sig.arg2;
    if (parameters.verbose) printf("Fetched adsorbate FF parameters\n");

    //
    // Construct array of framework particles, framework_atoms, for speed in energy computations
    //
    Particle_f * framework_atoms = (Particle_f *) malloc(framework.noatoms * sizeof(Particle_f));
    load_fast_particle_f_array(framework_atoms, framework, forcefield, parameters.epsilon_guest, parameters.sigma_guest);
    if (parameters.verbose) printf("Initialized framework_atoms array in host\n");

    parameters.num_blocks = NUM_BLOCKS; parameters.num_threads = NUM_THREADS;
	parameters.numinsertionsperthread = parameters.numinsertions / (NUM_BLOCKS * NUM_THREADS) + 1;
    if (parameters.verbose) printf("# blocks = %d, # threads = %d\n", parameters.num_blocks, parameters.num_threads);

    //
    // Write settings to output file
    //
    FILE * outputfile;
    char outputfilename[512];
    sprintf(outputfilename, "output_files/%s_%s_henry.out", parameters.frameworkname.c_str(), parameters.adsorbate.c_str());
    outputfile = fopen(outputfilename, "w");
    write_settings_to_outputfile(outputfile, parameters, framework, forcefield, framework_atoms);
    if (parameters.verbose) printf("Wrote info to outputfile\n");
  
	//
	// Set up random number generator on device
	//
    curandStateMtgp32 *devMTGPStates;
    mtgp32_kernel_params *devKernelParams;
      
    // Allocate space for prng states on device 
    CUDA_CALL(cudaMalloc((void **)&devMTGPStates, NUM_BLOCKS * sizeof(curandStateMtgp32)));

    // Allocate space for MTGP kernel parameters 
    CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));

    // Reformat from predefined parameter sets to kernel format, 
    // and copy kernel parameters to device memory               
    CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));

    // Initialize one state per thread block */
    CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, 
            mtgp32dc_params_fast_11213, devKernelParams, NUM_BLOCKS, 1234));

	//
    //  Move data to GPU device
    //
	Henry_stats * d_statistics; // collect statistics from device block
    CUDA_CALL(cudaMalloc((void **) & d_statistics, sizeof(Henry_stats) * NUM_BLOCKS)); 
	Particle_f * d_framework_atoms;
	CUDA_CALL(cudaMalloc((void **) & d_framework_atoms, framework.noatoms * sizeof(Particle_f)));
	CUDA_CALL(cudaMemcpy(d_framework_atoms, framework_atoms, framework.noatoms * sizeof(Particle_f), cudaMemcpyHostToDevice));
//	fprintf(outputfile, "    Size of framework atoms array = %f MB\n", framework.noatoms * sizeof(Particle_f) / (1024.0 * 1024.0));
	
    //
	//  RUN WIDOM INSERTIONS
	//
	double start_of_sim_time=read_timer();
	doHenry <<< NUM_BLOCKS, NUM_THREADS>>> (devMTGPStates, parameters, d_statistics, d_framework_atoms);
 
	cudaDeviceSynchronize();
	
    // get statistics
	Henry_stats * h_statistics= (Henry_stats *) malloc(NUM_BLOCKS * sizeof(Henry_stats));
    CUDA_CALL(cudaMemcpy(h_statistics, d_statistics, sizeof(Henry_stats) * NUM_BLOCKS,cudaMemcpyDeviceToHost));
 	
    double sim_time = read_timer() - start_of_sim_time;
   
    fprintf(outputfile, "\nHenry simulation results\n");
    fprintf(outputfile, "    Simulation time: %f s\n", sim_time);
	cudaDeviceSynchronize();
    if (parameters.verbose) printf("Copied data back from device to host\n");

	double henry_coeff = 0.0;
	double canonical_sum = 0.0;
	double weighted_energy_sum = 0.0;

	for (int i = 0; i < NUM_BLOCKS; i ++) {
		canonical_sum += h_statistics[i].canonical_sum;
		weighted_energy_sum += h_statistics[i].weighted_energy_sum;
		if (parameters.verbose) 
            printf("Block %d: canonical sum = %f\n", i, h_statistics[i].canonical_sum);
	}
	
	double ensemble_average_energy = weighted_energy_sum / canonical_sum;
	henry_coeff = canonical_sum / (NUM_THREADS * NUM_BLOCKS * parameters.numinsertionsperthread) / parameters.T / 8.314;
	fprintf(outputfile,  "    Temperature (K): %f\n", parameters.T);
	fprintf(outputfile,  "    <energy> (kJ/mol): %f\n", ensemble_average_energy * 8.314 /1000.0);
	fprintf(outputfile,  "    <energy> (K): %f\n", ensemble_average_energy);
	fprintf(outputfile,  "    Henry coefficient (mol/(kg-Pa)): %f\n", henry_coeff / framework.density);
	fprintf(outputfile,  "    Henry coefficient (mol/(m3-Pa)): %f\n", henry_coeff);
}
