/*
 * Writes Potential energy grid of adsorbate inside unit cell of nanoporous material
 */
#include <stdio.h>
#include <stdlib.h>
#include<cuda.h>
#include <sys/time.h>
#include <cuda_runtime.h>
using namespace std;
#include "datatypes.h"
#include "readsettings.h"
#include "Framework.h"
#include "Forcefield.h"
#include "write_settings_to_outputfile.h"
#include "computegridsheet.h"
#include "load_fast_particle_f_array.h"


#define SQRT_N_THREADS 16 // a block may have a max of 512 threads... so 16x16 is max.

// functions to ensure communication with GPU works ?
#define CUDA_CALL(x) do { cudaError_t error = x;  \
  if (error != cudaSuccess) {  \
  printf("Error at %s:%d - %s \n",__FILE__,__LINE__, cudaGetErrorString(error)); \
  return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
            return EXIT_FAILURE;}} while(0)

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
        printf("Run as:\n./writegrid structure_name AdsorbateID\n");
        exit(EXIT_FAILURE);
    }
    //
    //  Import settings
    //
    GridParameters parameters;
    parameters.frameworkname = argv[1];
    parameters.adsorbate = argv[2];
    parameters.adsorbateMW = GetAdsorbateMW(parameters.adsorbate);

    ReadSimulationInputFile(parameters);
    if (parameters.verbose) printf("Read simulation.input\n");
    // only need UC to be once the cutoff
    TripleInt uc_reps = ReadUnitCellReplicationFile(parameters.frameworkname, "once");
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
    parameters.N_framework_atoms = framework.noatoms;
    if (parameters.verbose) printf("Constructed Framework object\n");

    // grab sigma/epsilon of adsorbate
    PairDouble eps_sig = GrabGuestForceFieldParams(forcefield, parameters.adsorbate);
    parameters.epsilon_guest = eps_sig.arg1;
    parameters.sigma_guest = eps_sig.arg2;
    if (parameters.verbose) printf("Fetched adsorbate FF parameters\n");

    //
    // Construct array of framework particles, framework_atoms, for speed in energy computations
    //
    FrameworkParticle * framework_atoms = (FrameworkParticle *) malloc(framework.noatoms * sizeof(FrameworkParticle));
    LoadFastFrameworkParticleArray(framework_atoms, framework, forcefield, parameters.epsilon_guest, parameters.sigma_guest);
    if (parameters.verbose) printf("Initialized framework_atoms array in host\n");

    //
    // Construct grid
    //
    int N_x = static_cast<int>(ceil(framework.a / parameters.grid_resolution)); // size of grid
    int N_y = static_cast<int>(ceil(framework.b / parameters.grid_resolution));
    int N_z = static_cast<int>(ceil(framework.c / parameters.grid_resolution));
    parameters.N_x = N_x; parameters.N_y = N_y; parameters.N_z = N_z;

    // pointer array of fractional grid points
    double * x_f_gridpoints = (double *) malloc(N_x * sizeof(double));
    double * y_f_gridpoints = (double *) malloc(N_y * sizeof(double));
    double * z_f_gridpoints = (double *) malloc(N_z * sizeof(double));

    // fractional coordinate for a unit cell
    for (int i = 0; i < N_x; i++)
        x_f_gridpoints[i] = 1.0 * i / (N_x - 1);
    for (int i = 0; i < N_y; i++)
        y_f_gridpoints[i] = 1.0 * i / (N_y - 1);
    for (int i = 0; i < N_z; i++)
        z_f_gridpoints[i] = 1.0 * i / (N_z - 1);

    //
    //  Write settings to outputfile
    //
    FILE * outputfile;
    char outputfilename[512];
    sprintf(outputfilename, "output_files/%s_%s_grid.out", parameters.frameworkname.c_str(), parameters.adsorbate.c_str());
    outputfile = fopen(outputfilename, "w");
    WriteSettingsToOutputfile(outputfile, parameters, framework, forcefield, framework_atoms);
    if (parameters.verbose) printf("Wrote info to outputfile\n");

    //
    // PREPARE GRID FILE
    //
    FILE * gridfile;
    char gridfilename[256] = "data/grids/";
    strcat(gridfilename, framework.name.c_str());
    strcat(gridfilename,"_"); strcat(gridfilename, parameters.adsorbate.c_str());
    strcat(gridfilename,"_"); strcat(gridfilename, forcefield.name.c_str());
    if (parameters.gridoutputformat == "txt") {  // format I  made up for Henry coefficient and GCMC calcs
        gridfile = fopen(strcat(gridfilename, ".txt"), "w");
        fprintf(gridfile, "%d %d %d  = (N_x,N_y,N_z) grid points (grid is in fractional coords). Endpoints included.\n", N_x, N_y, N_z);
    }
    else if (parameters.gridoutputformat == "cube") // for visualization with VisIt
    {
        gridfile = fopen(strcat(gridfilename, ".cube"), "w");
        fprintf(gridfile, "\nThis is a grid file.\n");
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              framework.noatoms, 0.0, 0.0, 0.0); // give number of atoms
        // give little vectors that form a volume element
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              N_x, framework.t_matrix[0][0] / (N_x - 1), 0.0, 0.0);
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              N_y, framework.t_matrix[0][1]/ (N_y - 1), framework.t_matrix[1][1] / (N_y - 1), 0.0);
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              N_z, framework.t_matrix[0][2] / (N_z - 1), framework.t_matrix[1][2] / (N_z - 1), framework.t_matrix[2][2] / (N_z - 1));

        // write atoms to grid file
        double x, y, z;
      for (int i = 0; i < framework.noatoms; i++) {
          double atomic_mass = 1.0;
          int atomic_number = 1;
          // TODO get atomic numbers and atomic masses, for now just say they are all H ha ha
          HostFractionalToCartesian(framework.t_matrix, framework_atoms[i].x_f, framework_atoms[i].y_f, framework_atoms[i].z_f , x, y, z);
          fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf % 13.6lf\n", atomic_number, atomic_mass, x, y, z);
      }
      fprintf(gridfile," 1    1\n");
    }
    else {
        printf("Grid output format must be txt or cube\n");
        exit(EXIT_FAILURE);
    }
    if (parameters.verbose) printf("Initialized grid file\n");

    //
    // Parallelization strategy: pass sheets of the grid to the GPU at a time, sheets are defind by x = constant
    //
    // energies at zy grid sheet. entry k+j*N_z is the energy at point z=k*dz, y = j*dy
    double * h_zy_energies = (double *) malloc(N_z * N_y * sizeof(double));

    //
    //  Move data to GPU device; "d_" indicates this is data for the device
    //
    // Initialize memory for zy_energies on device, to be called and stored bck to zy_energies later
    double * d_zy_energies;
    CUDA_CALL(cudaMalloc((void **) & d_zy_energies, N_z * N_y * sizeof(double)));
    // Copy framework_atoms to device. All blocks share this.
    FrameworkParticle * d_framework_atoms;
    CUDA_CALL(cudaMalloc((void **) & d_framework_atoms, framework.noatoms * sizeof(FrameworkParticle)));
    CUDA_CALL(cudaMemcpy(d_framework_atoms, framework_atoms, framework.noatoms * sizeof(FrameworkParticle), cudaMemcpyHostToDevice));
    fprintf(outputfile, "    Size of framework atoms array = %f MB\n", framework.noatoms * sizeof(FrameworkParticle) / (1024.0 * 1024.0));
    // copy z_f and y_f grid points to device. The parallelization strategy is to pass sheets of x = constant, so this is not needed on the device.
    double * d_z_f_gridpoints;
    double * d_y_f_gridpoints;
    CUDA_CALL(cudaMalloc((void **) & d_z_f_gridpoints, N_z * sizeof(double)));
    CUDA_CALL(cudaMalloc((void **) & d_y_f_gridpoints, N_y * sizeof(double)));
    CUDA_CALL(cudaMemcpy(d_z_f_gridpoints, z_f_gridpoints, N_z * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_y_f_gridpoints, y_f_gridpoints, N_y * sizeof(double), cudaMemcpyHostToDevice));
    fprintf(outputfile, "    Size of grid sheet = %f MB\n", N_z * N_y * sizeof(double) / (1024.0 * 1024.0));
    if (parameters.verbose) printf("Copied framework_atoms, z_f/y_f grid points, and allocated zy_energies to GPU device\n");

    //
    // Write the grid
    //
    fprintf(outputfile, "    A block is %d by %d threads.\n", SQRT_N_THREADS, SQRT_N_THREADS);
    dim3 dimBlock(SQRT_N_THREADS, SQRT_N_THREADS); // size of block. making 2D thread block
    dim3 dimGrid(N_z / SQRT_N_THREADS + 1, N_y / SQRT_N_THREADS + 1);
    double t0 = ReadTimer();
    if (parameters.verbose) printf("Starting loop to write grid...\n# x-grid points: %d\n", N_x);

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++){
            parameters.t_matrix[i][j] = framework.t_matrix[i][j];
        }
    }// TODO: remove and use framework.t_matrix instead. right now it cant pass to cuda kernel without memory error...

    int count_grid_pts = 0;
    for (int i = 0; i < N_x; i++) {
        ComputeGridSheet <<<dimGrid, dimBlock>>> (d_z_f_gridpoints,
                d_y_f_gridpoints,
                d_zy_energies,
                d_framework_atoms,
                parameters,
                x_f_gridpoints[i]);
        CUDA_CALL( cudaPeekAtLastError() );
        CUDA_CALL( cudaDeviceSynchronize() );
        cudaDeviceSynchronize();

        // get energies from device
        CUDA_CALL(cudaMemcpy(h_zy_energies, d_zy_energies, N_z * N_y * sizeof(double) , cudaMemcpyDeviceToHost));
        cudaDeviceSynchronize();

        // write energies to file
        if (parameters.gridoutputformat=="cube") {
            for (int j=0; j < N_y; j++) {
                int count = 0;
                for(int k = 0; k < N_z; k++) {
                    fprintf(gridfile, "% 13.6E ", h_zy_energies[k + j * N_z] * 8.314 / 1000); // kJ/mol
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

        if (parameters.gridoutputformat=="txt") { // format I made up ^.^ TODO more efficient format?
            for (int j = 0; j < N_y; j++) {
                for(int k = 0; k < N_z; k++) {
                    count_grid_pts += 1;
                    fprintf(gridfile, "% 13.6E ", h_zy_energies[k + j * N_z]);
                    if ( k == (N_z - 1))
                        fprintf(gridfile, "\n"); // new line for every pencil of z's
                }
            }
        }
        if (parameters.verbose) printf("   Sheet %d out of %d completed.\n", i, N_x);
    } // end x loop

    double sim_time = ReadTimer() - t0;
    fprintf(outputfile, "    Time to write grid: %f s\n", sim_time);
    if (parameters.verbose) printf("Completed grid writing! Freeing up memory in GPU...\n");

    //
    // Free memory, close files
    //
    cudaFree(d_framework_atoms);
    cudaFree(d_zy_energies);
    cudaFree(d_z_f_gridpoints);
    cudaFree(d_y_f_gridpoints);
    free(framework_atoms);
    free(x_f_gridpoints); free(y_f_gridpoints); free(z_f_gridpoints);
    fclose(outputfile); fclose(gridfile);

}

