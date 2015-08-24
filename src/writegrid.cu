/*
 * Writes Potential energy grid of adsorbate inside unit cell of nanoporous material
 */
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
#include "computegridsheet.h"
#include "load_fast_particle_f_array.h"
#include "GetEwaldParams.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>


#define SQRT_N_THREADS 16 // a block may have a max of 512 threads... so 16x16 is max.

// functions to ensure communication with GPU works ?
#define CUDA_CALL(x) do { cudaError_t error = x;  \
  if (error != cudaSuccess) {  \
  printf("Error at %s:%d - %s \n",__FILE__,__LINE__, cudaGetErrorString(error)); \
  return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
            return EXIT_FAILURE;}} while(0)

thrust::host_vector<double> cross_product(thrust::host_vector<double> & a, thrust::host_vector<double> & b) {
    // return cross product of vectors a and b (3D only)
    thrust::host_vector<double> result(3);
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
}

double dot_product(thrust::host_vector<double> & a, thrust::host_vector<double> & b) {
    // return cross product of vectors a and b (3D only)
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double ReadTimer() {
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

void HostFractionalToCartesian(double t_matrix[][3],
        double x_f, double y_f,double z_f,
        double & x, double & y, double & z) {
    // compute Cartesian coordinates from fractional
    x = t_matrix[0][0] * x_f + t_matrix[0][1] * y_f + t_matrix[0][2] * z_f;
    y = t_matrix[1][0] * x_f + t_matrix[1][1] * y_f + t_matrix[1][2] * z_f;
    z = t_matrix[2][0] * x_f + t_matrix[2][1] * y_f + t_matrix[2][2] * z_f;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Run as:\n./writegrid structure_name adsorbate_bead\n");
        exit(EXIT_FAILURE);
    }

    //
    //  Import settings
    //
    GridParameters parameters;
    parameters.frameworkname = argv[1];
    parameters.adsorbatebead = argv[2];
    // Coulomb grid if adsorbate is Coulomb
    if (parameters.adsorbatebead == "Coulomb")
        parameters.Coulomb_grid_flag = true;
    else
        parameters.Coulomb_grid_flag = false;
//    parameters.adsorbatebeadMW = GetAdsorbateMW(parameters.adsorbatebead);

    ReadSimulationInputFile(parameters);
    if (parameters.verbose) 
        printf("Read simulation.input\n");
    if ((parameters.verbose) & (parameters.Coulomb_grid_flag))
        printf("----------------   Coloumb Grid. -----------------------\n");

    // Get unit cell replication factors.
    //    only need UC to be once the cutoff
    std::vector<int> uc_reps = ReadUnitCellReplicationFile(parameters.frameworkname, "once");
    parameters.replication_factor_a = uc_reps[0];
    parameters.replication_factor_b = uc_reps[1];
    parameters.replication_factor_c = uc_reps[2];
    if (parameters.verbose) 
        printf("Read .uc replication file\n");

    //
    // Construct forcefield and framework objects
    //
    Forcefield forcefield(parameters.forcefieldname);
    if (parameters.verbose) 
        printf("Constructed Forcefield object\n");
    
    Framework framework(parameters.frameworkname);
    parameters.N_framework_atoms = framework.noatoms;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++){
            parameters.t_matrix[i][j] = framework.t_matrix[i][j];
        }
    }
    parameters.volume_unitcell = framework.volume_unitcell;

    if (parameters.verbose) 
        printf("Constructed Framework object\n");

    // grab sigma/epsilon of adsorbate bead
    std::vector<double> eps_sig = GrabGuestForceFieldParams(forcefield, parameters.adsorbatebead);
    parameters.epsilon_guest = eps_sig[0];
    parameters.sigma_guest = eps_sig[1];
    if (parameters.verbose) 
        printf("Fetched adsorbate FF parameters\n");

    //
    // Construct array of framework particles, framework_atoms, for speed in energy computations
    //
    FrameworkParticle * framework_atoms = (FrameworkParticle *) malloc(framework.noatoms * sizeof(FrameworkParticle));
    LoadFastFrameworkParticleArray(framework_atoms, framework, forcefield, parameters.epsilon_guest, parameters.sigma_guest);
    if (parameters.verbose) 
        printf("Initialized framework_atoms array in host\n");
    
    //
    // Construct reciprocal lattice vectors on device for EWald summations
    //
    std::vector<int> all_ones(3); all_ones[0] = 1; all_ones[1] = 1; all_ones[2] = 1;
    EWaldParameters ew_params = GetEwaldParams(framework, all_ones, sqrt(parameters.r_cutoff_squared), 
                                                parameters.EWald_precision, parameters.verbose);
    thrust::device_vector<double> b1(3);
    thrust::device_vector<double> b2(3);
    thrust::device_vector<double> b3(3);
    // transfer from Ewald params to here
    for (int i = 0; i < 3; i++) {
        b1[i] = ew_params.b1[i];
        b2[i] = ew_params.b2[i];
        b3[i] = ew_params.b3[i];
    }

    //
    // Construct grid
    //
    parameters.N_x = static_cast<int>(ceil(framework.a / parameters.grid_resolution)); // size of grid
    parameters.N_y = static_cast<int>(ceil(framework.b / parameters.grid_resolution));
    parameters.N_z = static_cast<int>(ceil(framework.c / parameters.grid_resolution));

    // pointer array of fractional grid points
    double * x_f_gridpoints = (double *) malloc(parameters.N_x * sizeof(double));
    double * y_f_gridpoints = (double *) malloc(parameters.N_y * sizeof(double));
    double * z_f_gridpoints = (double *) malloc(parameters.N_z * sizeof(double));

    // fractional coordinate for a unit cell
    for (int i = 0; i < parameters.N_x; i++)
        x_f_gridpoints[i] = 1.0 * i / (parameters.N_x - 1);
    for (int i = 0; i < parameters.N_y; i++)
        y_f_gridpoints[i] = 1.0 * i / (parameters.N_y - 1);
    for (int i = 0; i < parameters.N_z; i++)
        z_f_gridpoints[i] = 1.0 * i / (parameters.N_z - 1);

    //
    //  Write settings to outputfile
    //
    FILE * outputfile;
    char outputfilename[512];
    sprintf(outputfilename, "output_files/%s_%s_grid.out", parameters.frameworkname.c_str(), parameters.adsorbatebead.c_str());
    outputfile = fopen(outputfilename, "w");
    WriteSettingsToOutputfile(outputfile, parameters, ew_params, framework, forcefield, framework_atoms);
    if (parameters.verbose) 
        printf("Wrote info to outputfile\n");

    //
    // PREPARE GRID FILES
    //
    FILE * gridfile;
    char gridfilename[512];
    if (parameters.gridoutputformat == "txt") {  // format I  made up for Henry coefficient and GCMC calcs
        if (! parameters.Coulomb_grid_flag)
            sprintf(gridfilename, "/home/corymsimon/sim_data/grids/vdW/%s_%s_%s.txt", framework.name.c_str(), parameters.adsorbatebead.c_str(), forcefield.name.c_str());
        else
            sprintf(gridfilename, "/home/corymsimon/sim_data/grids/Coulomb/%s.txt", framework.name.c_str());
        gridfile = fopen(gridfilename, "w");
        fprintf(gridfile, "%d %d %d  = (parameters.N_x,parameters.N_y,parameters.N_z)"
                                           "grid points (grid is in fractional coords)." 
                                           "Endpoints included.\n", parameters.N_x, parameters.N_y, parameters.N_z);
    }
    else if (parameters.gridoutputformat == "cube") {  // for visualization with VisIt
        if (! parameters.Coulomb_grid_flag)
            sprintf(gridfilename, "/sim_data/data/grids/vdW/%s_%s_%s.cube", framework.name.c_str(), parameters.adsorbatebead.c_str(), forcefield.name.c_str());
        else
            sprintf(gridfilename, "/sim_data/data/grids/Coulomb/%s.cube", framework.name.c_str());

        gridfile = fopen(gridfilename, "w");
        
        fprintf(gridfile, "\nThis is a grid file.\n");
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              0, 0.0, 0.0, 0.0); // give number of atoms
        // give little vectors that form a volume element
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              parameters.N_x, framework.t_matrix[0][0] / (parameters.N_x - 1), 0.0, 0.0);
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              parameters.N_y, framework.t_matrix[0][1] / (parameters.N_y - 1), framework.t_matrix[1][1] / (parameters.N_y - 1), 0.0);
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              parameters.N_z, framework.t_matrix[0][2] / (parameters.N_z - 1), framework.t_matrix[1][2] / (parameters.N_z - 1), framework.t_matrix[2][2] / (parameters.N_z - 1));
    }
    else {
        printf("Grid output format must be txt or cube\n");
        exit(EXIT_FAILURE);
    }
    if (parameters.verbose) printf("Initialized grid file\n");
    
    //
    // Parallelization strategy: pass sheets of the grid to the GPU at a time, sheets are defind by x = constant
    //
    // energies at zy grid sheet. entry k+j*parameters.N_z is the energy at point z=k*dz, y = j*dy
    double * h_zy_energies = (double *) malloc(parameters.N_z * parameters.N_y * sizeof(double));

    //
    //  Move data to GPU device; "d_" indicates this is data for the device
    //
    // Initialize memory for zy_energies on device, to be called and stored bck to zy_energies later
    double * d_zy_energies;  // van der waals energies
    CUDA_CALL(cudaMalloc((void **) & d_zy_energies, parameters.N_z * parameters.N_y * sizeof(double)));

    // Copy framework_atoms to device. All blocks share this.
    FrameworkParticle * d_framework_atoms;
    CUDA_CALL(cudaMalloc((void **) & d_framework_atoms, framework.noatoms * sizeof(FrameworkParticle)));
    CUDA_CALL(cudaMemcpy(d_framework_atoms, framework_atoms, framework.noatoms * sizeof(FrameworkParticle), cudaMemcpyHostToDevice));
    fprintf(outputfile, "    Size of framework atoms array = %f MB\n", framework.noatoms * sizeof(FrameworkParticle) / (1024.0 * 1024.0));
    // copy z_f and y_f grid points to device. The parallelization strategy is to pass sheets of x = constant, so this is not needed on the device.
    double * d_z_f_gridpoints;
    double * d_y_f_gridpoints;
    CUDA_CALL(cudaMalloc((void **) & d_z_f_gridpoints, parameters.N_z * sizeof(double)));
    CUDA_CALL(cudaMalloc((void **) & d_y_f_gridpoints, parameters.N_y * sizeof(double)));
    CUDA_CALL(cudaMemcpy(d_z_f_gridpoints, z_f_gridpoints, parameters.N_z * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_y_f_gridpoints, y_f_gridpoints, parameters.N_y * sizeof(double), cudaMemcpyHostToDevice));
    fprintf(outputfile, "    Size of grid sheet = %f MB\n", parameters.N_z * parameters.N_y * sizeof(double) / (1024.0 * 1024.0));
    if (parameters.verbose) 
        printf("Copied framework_atoms, z_f/y_f grid points, and allocated zy_energies to GPU device\n");

    //
    // Write the grid
    //
    fprintf(outputfile, "    A block is %d by %d threads.\n", SQRT_N_THREADS, SQRT_N_THREADS);
    dim3 dimBlock(SQRT_N_THREADS, SQRT_N_THREADS); // size of block. making 2D thread block
    dim3 dimGrid(parameters.N_z / SQRT_N_THREADS + 1, parameters.N_y / SQRT_N_THREADS + 1);
    double t0 = ReadTimer();
    if (parameters.verbose) printf("Starting loop to write grid...\n# x-grid points: %d\n", parameters.N_x);

    int count_grid_pts = 0;

    for (int i = 0; i < parameters.N_x; i++) {
        if (parameters.verbose) 
            printf("   Sheet %d out of %d...\n", i, parameters.N_x);
        // van der Waals
        if ((! parameters.Coulomb_grid_flag)) {
            ComputevdWGridSheet <<<dimGrid, dimBlock>>> (d_z_f_gridpoints,
                    d_y_f_gridpoints,
                    d_zy_energies,
                    d_framework_atoms,
                    parameters,
                    x_f_gridpoints[i]);
            CUDA_CALL( cudaPeekAtLastError() );
            CUDA_CALL( cudaDeviceSynchronize() );
            
            // get energies from device
            CUDA_CALL(cudaMemcpy(h_zy_energies, d_zy_energies, parameters.N_z * parameters.N_y * sizeof(double) , cudaMemcpyDeviceToHost));
            cudaDeviceSynchronize();
        }
        
        // Coulomb
        if (parameters.Coulomb_grid_flag) {
            ComputeCoulombGridSheet <<<dimGrid, dimBlock>>> (d_z_f_gridpoints,
                    d_y_f_gridpoints,
                    d_zy_energies,
                    d_framework_atoms,
                    parameters,
                    ew_params,
                    thrust::raw_pointer_cast( &b1[0]),
                    thrust::raw_pointer_cast( &b2[0]),
                    thrust::raw_pointer_cast( &b3[0]),
                    x_f_gridpoints[i]);
            CUDA_CALL( cudaPeekAtLastError() );
            CUDA_CALL( cudaDeviceSynchronize() );
        
            // get energies from device
            CUDA_CALL(cudaMemcpy(h_zy_energies, d_zy_energies, parameters.N_z * parameters.N_y * sizeof(double) , cudaMemcpyDeviceToHost));
            cudaDeviceSynchronize();
        }

        // write energies to file
        if (parameters.gridoutputformat=="cube") {
            // if cube file, write total energy in kJ/mol
            for (int j = 0; j < parameters.N_y; j++) {
                int count = 0;
                for(int k = 0; k < parameters.N_z; k++) {
                    double E_here = h_zy_energies[k + j * parameters.N_z]; // K
                    fprintf(gridfile, "% 13.6E ", E_here * 8.314 / 1000); // kJ/mol
                    count ++;
                    if (count == 6) {
                        fprintf(gridfile, "\n");
                        count = 0; // reset counter
                    }
                    count_grid_pts ++;
                }
                fprintf(gridfile, "\n"); //new line after z over
            }
        }

        if (parameters.gridoutputformat=="txt") { 
            // format I made up ^.^ TODO more efficient format?
            // separate grid for vdW and Coulomb
            for (int j = 0; j < parameters.N_y; j++) {
                for(int k = 0; k < parameters.N_z; k++) {
                    count_grid_pts += 1;
                    fprintf(gridfile, "% 13.6E ", h_zy_energies[k + j * parameters.N_z]);  // K
                    if (k == (parameters.N_z - 1)) {
                        // new line for every pencil of z's
                        fprintf(gridfile, "\n"); 
                    }
                }
            }
        }
    } // end x loop
    assert(count_grid_pts == (parameters.N_x * parameters.N_y * parameters.N_z));

    double sim_time = ReadTimer() - t0;
    fprintf(outputfile, "    Time to write grid: %f s\n", sim_time);
    if (parameters.verbose) 
        printf("Completed grid writing! Freeing up memory in GPU...\n");

    //
    // Free memory, close files
    //
    cudaFree(d_framework_atoms);
    
    cudaFree(d_zy_energies);
    
    cudaFree(d_z_f_gridpoints);
    cudaFree(d_y_f_gridpoints);
    
    free(framework_atoms);
    
    free(x_f_gridpoints); 
    free(y_f_gridpoints); 
    free(z_f_gridpoints);
    
    fclose(outputfile); 
    fclose(gridfile);
}
