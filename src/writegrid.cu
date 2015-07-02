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
    parameters.adsorbate = argv[2];
//    parameters.adsorbateMW = GetAdsorbateMW(parameters.adsorbate);

    ReadSimulationInputFile(parameters);
    if (parameters.verbose) 
        printf("Read simulation.input\n");
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
    if (parameters.verbose) 
        printf("Constructed Framework object\n");

    //
    // Construct reciprocal lattice vectors on device for EWald summations
    //
    // primitive lattice vectors a
    thrust::host_vector<double> a1(3);
    thrust::host_vector<double> a2(3);
    thrust::host_vector<double> a3(3);
    HostFractionalToCartesian(parameters.t_matrix,
        1.0, 0.0, 0.0,
        a1[0], a1[1], a1[2]);
    HostFractionalToCartesian(parameters.t_matrix,
        0.0, 1.0, 0.0,
        a2[0], a2[1], a2[2]);
    HostFractionalToCartesian(parameters.t_matrix,
        0.0, 0.0, 1.0,
        a3[0], a3[1], a3[2]);
    
    // reciprocal lattice vectors b
    thrust::device_vector<double> b1(3);
    thrust::device_vector<double> b2(3);
    thrust::device_vector<double> b3(3);
    // cross products
    thrust::host_vector<double> a2_a3 = cross_product(a2, a3);
    thrust::host_vector<double> a3_a1 = cross_product(a3, a1);
    thrust::host_vector<double> a1_a2 = cross_product(a1, a2);

    for (int i = 0; i < 3; i++) {
        // factor of 2pi?
        b1[i] = a2_a3[i] / dot_product(a1, a2_a3);
        b2[i] = a3_a1[i] / dot_product(a2, a3_a1);
        b3[i] = a1_a2[i] / dot_product(a3, a1_a2);
    }

    if (parameters.verbose) {
        printf("Primitive lattice vectors:\n");
        std::cout << "  a1 = (" << a1[0] << ", " << a1[1] << ", " << a1[2] << ")" << std::endl;
        std::cout << "  a2 = (" << a2[0] << ", " << a2[1] << ", " << a2[2] << ")" << std::endl;
        std::cout << "  a3 = (" << a3[0] << ", " << a3[1] << ", " << a3[2] << ")" << std::endl;
        printf("Reciprocal lattice vectors:\n");
        std::cout << "  b1 = (" << b1[0] << ", " << b1[1] << ", " << b1[2] << ")" << std::endl;
        std::cout << "  b2 = (" << b2[0] << ", " << b2[1] << ", " << b2[2] << ")" << std::endl;
        std::cout << "  b3 = (" << b3[0] << ", " << b3[1] << ", " << b3[2] << ")" << std::endl;
    }
        
    

    // grab sigma/epsilon of adsorbate bead
    std::vector<double> eps_sig = GrabGuestForceFieldParams(forcefield, parameters.adsorbate);
    parameters.epsilon_guest = eps_sig[0];
    parameters.sigma_guest = eps_sig[1];
    if (parameters.verbose) 
        printf("Fetched adsorbate FF parameters\n");

    //
    // Construct array of framework particles, framework_atoms, for speed in energy computations
    //
    FrameworkParticle * framework_atoms = (FrameworkParticle *) malloc(framework.noatoms * sizeof(FrameworkParticle));
    LoadFastFrameworkParticleArray(framework_atoms, framework, forcefield, parameters.epsilon_guest, parameters.sigma_guest);
    if (parameters.verbose) printf("Initialized framework_atoms array in host\n");

    //
    // Construct grid
    //
    int N_x, N_y, N_z;
    N_x = static_cast<int>(ceil(framework.a / parameters.grid_resolution)); // size of grid
    N_y = static_cast<int>(ceil(framework.b / parameters.grid_resolution));
    N_z = static_cast<int>(ceil(framework.c / parameters.grid_resolution));
    parameters.N_x = N_x; 
    parameters.N_y = N_y; 
    parameters.N_z = N_z;

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
    if (parameters.verbose) 
        printf("Wrote info to outputfile\n");

    //
    // PREPARE GRID FILE
    //
    FILE * gridfile;
    char gridfilename[512];
    if (parameters.gridoutputformat == "txt") {  // format I  made up for Henry coefficient and GCMC calcs
        sprintf(gridfilename, "data/grids/%s_%s_%s.txt", framework.name.c_str(), parameters.adsorbate.c_str(), forcefield.name.c_str());
        gridfile = fopen(gridfilename, "w");
        fprintf(gridfile, "%d %d %d  = (N_x,N_y,N_z) grid points (grid is in fractional coords). Endpoints included.\n", N_x, N_y, N_z);
    }
    else if (parameters.gridoutputformat == "cube") {  // for visualization with VisIt
        sprintf(gridfilename, "data/grids/%s_%s_%s.cube", framework.name.c_str(), parameters.adsorbate.c_str(), forcefield.name.c_str());

        gridfile = fopen(gridfilename, "w");
        
        fprintf(gridfile, "\nThis is a grid file.\n");
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              0, 0.0, 0.0, 0.0); // give number of atoms
        // give little vectors that form a volume element
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              N_x, framework.t_matrix[0][0] / (N_x - 1), 0.0, 0.0);
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              N_y, framework.t_matrix[0][1] / (N_y - 1), framework.t_matrix[1][1] / (N_y - 1), 0.0);
        fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n",
              N_z, framework.t_matrix[0][2] / (N_z - 1), framework.t_matrix[1][2] / (N_z - 1), framework.t_matrix[2][2] / (N_z - 1));
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
    double * h_zy_vdW_energies = (double *) malloc(N_z * N_y * sizeof(double));
    double * h_zy_Coulomb_energies = (double *) malloc(N_z * N_y * sizeof(double));

    //
    //  Move data to GPU device; "d_" indicates this is data for the device
    //
    // Initialize memory for zy_energies on device, to be called and stored bck to zy_energies later
    double * d_zy_vdW_energies;  // van der waals energies
    CUDA_CALL(cudaMalloc((void **) & d_zy_vdW_energies, N_z * N_y * sizeof(double)));
    double * d_zy_Coulomb_energies;  // coulomb energies
    CUDA_CALL(cudaMalloc((void **) & d_zy_Coulomb_energies, N_z * N_y * sizeof(double)));

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
    if (parameters.verbose) 
        printf("Copied framework_atoms, z_f/y_f grid points, and allocated zy_energies to GPU device\n");

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
//            printf("x_F=%f\n", x_f_gridpoints[i]);
        ComputeGridSheet <<<dimGrid, dimBlock>>> (d_z_f_gridpoints,
                d_y_f_gridpoints,
                d_zy_vdW_energies,
                d_zy_Coulomb_energies,
                d_framework_atoms,
                parameters,
                x_f_gridpoints[i]);
        CUDA_CALL( cudaPeekAtLastError() );
        CUDA_CALL( cudaDeviceSynchronize() );
        cudaDeviceSynchronize();

        // get energies from device
        CUDA_CALL(cudaMemcpy(h_zy_vdW_energies, d_zy_vdW_energies, N_z * N_y * sizeof(double) , cudaMemcpyDeviceToHost));
        cudaDeviceSynchronize();
        CUDA_CALL(cudaMemcpy(h_zy_Coulomb_energies, d_zy_Coulomb_energies, N_z * N_y * sizeof(double) , cudaMemcpyDeviceToHost));
        cudaDeviceSynchronize();
        
//            printf("\n\n Host:\n");
//            for (int kk = 0; kk < N_z; kk++) {
//                for (int jj=0; jj<N_y;jj++){
//                    printf("E[%d]=%f\n", kk+ jj*N_z, h_zy_energies[kk+ jj*N_z]);
//                }
//            }
//            exit(EXIT_FAILURE);

        // write energies to file
        if (parameters.gridoutputformat=="cube") {
            for (int j = 0; j < N_y; j++) {
                int count = 0;
                for(int k = 0; k < N_z; k++) {
                    fprintf(gridfile, "% 13.6E ", h_zy_vdW_energies[k + j * N_z] * 8.314 / 1000); // kJ/mol
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
                    fprintf(gridfile, "% 13.6E ", h_zy_vdW_energies[k + j * N_z]);
                    if ( k == (N_z - 1))
                        fprintf(gridfile, "\n"); // new line for every pencil of z's
                }
            }
        }
        if (parameters.verbose) printf("   Sheet %d out of %d completed.\n", i, N_x);
    } // end x loop
    assert(count_grid_pts == (N_x * N_y * N_z));

    double sim_time = ReadTimer() - t0;
    fprintf(outputfile, "    Time to write grid: %f s\n", sim_time);
    if (parameters.verbose) 
        printf("Completed grid writing! Freeing up memory in GPU...\n");

    //
    // Free memory, close files
    //
    cudaFree(d_framework_atoms);
    
    cudaFree(d_zy_vdW_energies);
    cudaFree(d_zy_Coulomb_energies);
    
    cudaFree(d_z_f_gridpoints);
    cudaFree(d_y_f_gridpoints);
    
    free(framework_atoms);
    
    free(x_f_gridpoints); 
    free(y_f_gridpoints); 
    free(z_f_gridpoints);
    
    fclose(outputfile); 
    fclose(gridfile);
}
