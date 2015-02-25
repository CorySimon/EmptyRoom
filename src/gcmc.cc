//
//  GCMC code (mu, V, T) simulations.
//
#include <chrono>
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include <limits>
#include <complex>
#include <cstring>
#include<math.h>
#include<random>
#include<assert.h>
#include<cstdlib> // for "exit"
#include<sstream> // string stream
#include<vector>
#include "datatypes.h"
#include "Framework.h"
#include "Forcefield.h"
#include <sys/time.h>
#include "readsettings.h"
#include "Framework.h"
#include "Forcefield.h"
#include "write_settings_to_outputfile.h"
#include <sys/time.h>
#include "pocketblocking.h"
#define min_r .0000000000000001  // don't want to divide by zero...
#define MAX_GUESTS 2000 // max number of guests in an array

void InitializeGCMCStats(GCMCStats & stats) {
    // initializes GCMC statistics to zero
    stats.N_move_trials = 0; stats.N_insertion_trials = 0; stats.N_deletion_trials = 0; stats.N_ID_swap_trials = 0;
    stats.N_insertions = 0; stats.N_deletions = 0; stats.N_moves = 0; stats.N_ID_swaps = 0;
    stats.N_samples = 0;
    stats.guest_guest_energy_avg = 0.0; stats.framework_guest_energy_avg= 0.0;
    stats.N_g_avg[0] = 0.0; stats.N_g_avg[1] = 0.0;
    stats.N_g2_avg[0] = 0.0; stats.N_g2_avg[1] = 0.0;
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

void FractionalToCartesian(double T_matrix[3][3],
    double x_f, double y_f, double z_f,
    double & x, double & y, double & z) 
{ 
    // compute cartesian coordinates from fractional
    x = T_matrix[0][0] * x_f + T_matrix[0][1] * y_f + T_matrix[0][2] * z_f;
    y = T_matrix[1][0] * x_f + T_matrix[1][1] * y_f + T_matrix[1][2] * z_f;
    z = T_matrix[2][0] * x_f + T_matrix[2][1] * y_f + T_matrix[2][2] * z_f;
}

void CartesianToFractional(double inv_T_matrix[3][3],
    double & x_f, double & y_f, double & z_f,
    double x, double y, double z)
{ 
    // cartesian to fractional
    x_f = inv_T_matrix[0][0] * x + inv_T_matrix[0][1] * y + inv_T_matrix[0][2] * z;
    y_f = inv_T_matrix[1][0] * x + inv_T_matrix[1][1] * y + inv_T_matrix[1][2] * z;
    z_f = inv_T_matrix[2][0] * x + inv_T_matrix[2][1] * y + inv_T_matrix[2][2] * z;
}

double WrapToInterval(double x, double z) {
    // for applying periodic bc's
    return x - z * floor(x / z);
}

//
// Compute potential energy of guest molecule id "which", the contribution from other guests.
//
double GuestGuestEnergy(int N_g, 
    int which, 
    GuestParticle * guests, 
    GCMCParameters parameters)
{
    double E_gg_ = 0.0; // initiate
    for (int k = 0 ; k < N_g; k++ ) { 
        // nearest particle image
        if (k == which) 
            continue; // do not include self interation orz
        // use nearest image
        double dx_f = guests[k].x_f - guests[which].x_f;
        double dy_f = guests[k].y_f - guests[which].y_f;
        double dz_f = guests[k].z_f - guests[which].z_f;
        dx_f = dx_f - parameters.replication_factor_a * round(dx_f / parameters.replication_factor_a);
        dy_f = dy_f - parameters.replication_factor_b * round(dy_f / parameters.replication_factor_b);
        dz_f = dz_f - parameters.replication_factor_c * round(dz_f / parameters.replication_factor_c);
//      assert (dx_f < 0.5 * parameters.replication_factor_a);
//      assert (dy_f < 0.5 * parameters.replication_factor_b);
//      assert (dz_f < 0.5 * parameters.replication_factor_c);
//      assert (dx_f > -0.5 * parameters.replication_factor_a);
//      assert (dy_f > -0.5 * parameters.replication_factor_b);
//      assert (dz_f > -0.5 * parameters.replication_factor_c);
        double dx, dy, dz;
        FractionalToCartesian(parameters.t_matrix,
            dx_f, dy_f, dz_f,
            dx, dy, dz);
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < min_r) 
            return 100000000000000.0; // overwrite if too small
        if (r2 < parameters.r_cutoff_squared) {
            double sigma_over_r_sixth = pow(parameters.sigma_squared_matrix[guests[which].type][guests[k].type] / r2, 3.0); // TODO: make sig2
            E_gg_ += 4.0 * parameters.epsilon_matrix[guests[which].type][guests[k].type] * sigma_over_r_sixth * (sigma_over_r_sixth - 1.0); // energy of framework atom id with guest_atom
        }
    }
    return E_gg_;
}

//
// Calculate total system guest-guest energy
// 
double GuestGuestTotalEnergy(int N_g, 
       GuestParticle * guests, 
       GCMCParameters parameters)
{
    double E_gg__ = 0.0;
    for (int i = 0; i < N_g; i ++) {
        E_gg__ += GuestGuestEnergy(N_g, i, guests, parameters); 
    }
    return E_gg__ / 2.0; // 2 is for double counting
}

//
//  Get index of flattened, 1D pointer array that corresponds to grid point index (i, j, k)
//
int GetEnergyGridIndex(int i, int j, int k, GridInfo grid_info) {
    // returns index in energy_grid corresponding to gridpoint index i,j,k
    return k + j * grid_info.N_z + i * grid_info.N_y * grid_info.N_z;
}

//
// Interpolate energy grid to get guest-framework energy
//  
double GuestFrameworkEnergy(double x_f_, double y_f_, double z_f_,
            GridInfo grid_info,
            double * energy_grid)
{
        // reflect fractional coords in [0,1] for grid
//      assert ((x_f_ > 0.0) & (y_f_ > 0.0) & (z_f_ > 0.0));
//      assert ((x_f_ < 1.0) & (y_f_ < 1.0) & (z_f_ < 1.0));
        //  FORMAT OF GRID: energy_grid[k+j*N_z+i*N_y*N_z]
        
        // define indices of 8 grid points, lower ones are:
        int i_x_low = floor(x_f_ / grid_info.dx_f);
        int i_y_low = floor(y_f_ / grid_info.dy_f);
        int i_z_low = floor(z_f_ / grid_info.dz_f);
        // trilinear interpolation http://en.wikipedia.org/wiki/Trilinear_interpolation
        
        // difference between our point and the vertices
        double x_d = (x_f_ - i_x_low * grid_info.dx_f) / grid_info.dx_f;
        double y_d = (y_f_ - i_y_low * grid_info.dy_f) / grid_info.dy_f;
        double z_d = (z_f_ - i_z_low * grid_info.dz_f) / grid_info.dz_f;
        
        // smash cube in x direction
        double c00 = energy_grid[GetEnergyGridIndex(i_x_low,i_y_low  ,i_z_low  ,grid_info)] * (1.0 - x_d) + x_d*energy_grid[GetEnergyGridIndex(i_x_low+1,i_y_low  ,i_z_low  ,grid_info)];
        double c10 = energy_grid[GetEnergyGridIndex(i_x_low,i_y_low+1,i_z_low  ,grid_info)] * (1.0 - x_d) + x_d*energy_grid[GetEnergyGridIndex(i_x_low+1,i_y_low+1,i_z_low  ,grid_info)];
        double c01 = energy_grid[GetEnergyGridIndex(i_x_low,i_y_low  ,i_z_low+1,grid_info)] * (1.0 - x_d) + x_d*energy_grid[GetEnergyGridIndex(i_x_low+1,i_y_low  ,i_z_low+1,grid_info)];
        double c11 = energy_grid[GetEnergyGridIndex(i_x_low,i_y_low+1,i_z_low+1,grid_info)] * (1.0 - x_d) + x_d*energy_grid[GetEnergyGridIndex(i_x_low+1,i_y_low+1,i_z_low+1,grid_info)];
        
        // further smash cube in y direction
        double c0 = c00 * (1.0 - y_d) + c10 * y_d;
        double c1 = c01 * (1.0 - y_d) + c11 * y_d;
        
        // finally, linear interpolation in z direction
        return c0 * (1 - z_d) + c1 * z_d;
}

//
// Calculate total system framework-guest energy
//
double GuestFrameworkTotalEnergy(int N_g, 
        double * energy_grid0, 
        double * energy_grid1, 
        GuestParticle * guests,
        GridInfo grid_info)
{
    double E_gf__ = 0.0;
    for (int i = 0; i < N_g; i ++) {
        if (guests[i].type == 0)
            E_gf__ += GuestFrameworkEnergy(WrapToInterval(guests[i].x_f,1.0), WrapToInterval(guests[i].y_f,1.0), WrapToInterval(guests[i].z_f,1.0),
                        grid_info, energy_grid0);
        if (guests[i].type == 1)
            E_gf__ += GuestFrameworkEnergy(WrapToInterval(guests[i].x_f,1.0), WrapToInterval(guests[i].y_f,1.0), WrapToInterval(guests[i].z_f,1.0),
                        grid_info, energy_grid1);
    }
    return E_gf__;
}

//
// Write guest positions to file
//
void WriteGuestPostionsToFile(FILE * positionfile, GuestParticle * guests, int N_g, GCMCParameters parameters) {
    for (int i = 0; i < N_g ; i++) {
        fprintf(positionfile, "%s %f %f %f\n", 
                parameters.adsorbate[guests[i].type].c_str(),
                guests[i].x, guests[i].y, guests[i].z);
    }
}

int main(int argc, char *argv[])
{
    if (! ((argc == 4) | (argc == 6))) {
        printf("Run as ./gcmc $structure $adsorbate0 $fugacity0(Pa) $adsorbate1 $fugactiy1(Pa)\nAdsorbate1 stuff is optional\n");
        exit(EXIT_FAILURE);
    }
    
    GCMCParameters parameters;
    // read arguments to get adsorbate and fugacity
    parameters.frameworkname = argv[1];
    parameters.adsorbate[0] = argv[2];
    parameters.fugacity[0] = atof(argv[3]);
    parameters.numadsorbates = 1; // overwite later if two
    if (argc == 6) {
        parameters.adsorbate[1] = argv[4];
        parameters.fugacity[1] = atof(argv[5]);
        parameters.numadsorbates = 2;
    }

    ReadSimulationInputFile(parameters);
    if (parameters.verbose) printf("Read simulation.input\n");
    if ((parameters.numadsorbates == 1) & (parameters.p_identity_change > 0.0)) {
        printf("Only 1 adsorbate and ID swap probability > 1... Make it zero.\n");
        exit(EXIT_FAILURE);
    }
    
    // uc needs to be at least twice the cutoff radius for only methane within r_c to be within cutoff
    TripleInt uc_reps = ReadUnitCellReplicationFile(parameters.frameworkname, "twice");
    parameters.replication_factor_a = uc_reps.arg1;
    parameters.replication_factor_b = uc_reps.arg2;
    parameters.replication_factor_c = uc_reps.arg3;
    if (parameters.verbose) printf("Read .uc replication file\n");
    
    parameters.adsorbateMW[0] = GetAdsorbateMW(parameters.adsorbate[0]);
    if (parameters.numadsorbates == 2)
        parameters.adsorbateMW[1] = GetAdsorbateMW(parameters.adsorbate[1]);

    //
    // Construct forcefield and framework objects
    //
    Forcefield forcefield(parameters.forcefieldname);
    if (parameters.verbose) printf("Constructed Forcefield object\n");
    Framework framework(parameters.frameworkname);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++){
            parameters.t_matrix[i][j] = framework.t_matrix[i][j];
            parameters.inv_t_matrix[i][j] = framework.inv_t_matrix[i][j];
        }
    }
    parameters.N_framework_atoms = framework.noatoms;
    if (parameters.verbose) printf("Constructed Framework object\n");
    
    //
    // Get adsorbate LJ params and store in epsilon and sigma_squared matrices
    //
    PairDouble eps_sig = GrabGuestForceFieldParams(forcefield, parameters.adsorbate[0]);
    parameters.epsilon_matrix[0][0] = eps_sig.arg1; 
    parameters.sigma_squared_matrix[0][0] = eps_sig.arg2 * eps_sig.arg2;
    if (parameters.numadsorbates == 2) {
        eps_sig = GrabGuestForceFieldParams(forcefield, parameters.adsorbate[1]);
        parameters.epsilon_matrix[1][1] = eps_sig.arg1; 
        parameters.sigma_squared_matrix[1][1] = eps_sig.arg2 * eps_sig.arg2;
        // mixing rules
        parameters.epsilon_matrix[0][1] = sqrt(parameters.epsilon_matrix[0][0] * parameters.epsilon_matrix[1][1]);
        parameters.epsilon_matrix[1][0] = parameters.epsilon_matrix[0][1]; // symmetric
        double mixed_sigma = (sqrt(parameters.sigma_squared_matrix[0][0]) + sqrt(parameters.sigma_squared_matrix[1][1])) / 2.0;
        parameters.sigma_squared_matrix[0][1] = mixed_sigma * mixed_sigma;
        parameters.sigma_squared_matrix[1][0] = mixed_sigma * mixed_sigma;
    }
    if (parameters.verbose) printf("Fetched adsorbate FF parameters\n");

    //
    // Import energy grid
    // + pocket blocking if enabled.
    //
    GridInfo grid_info; // for storing grid info
    double * energy_grid0; double * energy_grid1;
    bool pocket_block_verbose_mode = false; // writes grid before and after
    for (int n_c = 0; n_c < parameters.numadsorbates; n_c ++) {
        int N_x_temp, N_y_temp, N_z_temp;
        if (parameters.verbose) printf("Importing energy grid %d\n", n_c);
        char gridfilename[512];
        sprintf(gridfilename, "data/grids/%s_%s_%s.txt", framework.name.c_str(), parameters.adsorbate[n_c].c_str(), forcefield.name.c_str());

        std::ifstream gridfile(gridfilename); // input file stream
        if (gridfile.fail()) {
            printf("Grid file %s could not be loaded...\n", gridfilename);
            exit(EXIT_FAILURE);
        }
        std::string line;
        getline(gridfile, line);
        std::istringstream this_line(line);
        this_line >> N_x_temp >> N_y_temp >> N_z_temp;
       
        if (n_c == 0) {
            // load grid_info
            grid_info.N_x = N_x_temp; grid_info.N_y = N_y_temp; grid_info.N_z = N_z_temp;
            grid_info.numtotalpts = N_x_temp * N_y_temp * N_z_temp;
        }
        if (n_c == 1) {
            // assert 2nd energy grid has same resolution...
            assert(N_x_temp = grid_info.N_x);
            assert(N_y_temp = grid_info.N_y);
            assert(N_z_temp = grid_info.N_z);
        }

        // create grid import format
        if (n_c == 0)
            energy_grid0 = (double *) calloc(grid_info.numtotalpts, sizeof(double));
        if (n_c == 1)
            energy_grid1 = (double *) calloc(grid_info.numtotalpts, sizeof(double));

        int i = 0; int j = 0;
        while(getline(gridfile,line)) {
            std::istringstream this_line(line);
            for (int k = 0; k < grid_info.N_z; k ++) {  // each line is a pencil of z's
                int index_here = k + j * grid_info.N_z + i * grid_info.N_y * grid_info.N_z;
                if (n_c == 0)
                    this_line >> energy_grid0[index_here];
                if (n_c == 1)
                    this_line >> energy_grid1[index_here];
            }
            j += 1; // each line is a pencil of z's for a particular y... so update y index
            if (j == grid_info.N_y) {
                i += 1; // at the end of the z-y sheet, we go to a new x.
                j = 0; // j goes back to zero
            }
        }
        assert(i == grid_info.N_x); // assert that we are at the end of the file
        assert(j == 0);
        gridfile.close();
        // grid spacing in fractional coordinates  (how it was done in writegrid)
        grid_info.dx_f = 1.0/(grid_info.N_x - 1); // grid spacings
        grid_info.dy_f = 1.0/(grid_info.N_y - 1);
        grid_info.dz_f = 1.0/(grid_info.N_z - 1);
        if (parameters.verbose) printf("energy grid %d imported successfully.\n",n_c);
        
        //
        // Flood fill/ pocket blocking, if enabled
        //
        if (parameters.pocketblocking) {
            double time_before = ReadTimer();
            if (pocket_block_verbose_mode) {
                printf("Pocket blocking beginning. Write a cube for before\n");
                if (n_c == 0)
                    WriteCube("before_blocking_0", framework, parameters, energy_grid0, grid_info); //just so you can see the grid before ...
                if (n_c == 1)
                    WriteCube("before_blocking_1", framework, parameters, energy_grid1, grid_info); //just so you can see the grid before ...
            }
            if (n_c == 0)
                grid_info.numpockets[0] = FindAndBlockPockets(energy_grid0, grid_info, parameters.T, parameters);
            if (n_c == 1)
                grid_info.numpockets[1] = FindAndBlockPockets(energy_grid1, grid_info, parameters.T, parameters);
            // energy_grid is passed as a poitner, so it is modified now if accessible pockets were there.
            if (pocket_block_verbose_mode) {
                printf("Pocket blocking finished. Write a cube for after\n");
                double time_after = ReadTimer();
                printf("Time spent to find and block pockets: %f s\n", time_after - time_before);
                if (n_c == 0)
                    WriteCube("after_blocking_0", framework, parameters, energy_grid0, grid_info); // ... and after pocket blocking
                if (n_c == 1)
                    WriteCube("after_blocking_1", framework, parameters, energy_grid1, grid_info); // ... and after pocket blocking
            }
        }
    }
    
    //
    // Initialize stats and guests array (includes both components via "type" attribute)
    //
    GCMCStats stats;
    InitializeGCMCStats(stats);
    
    double volume = framework.volume_unitcell * parameters.replication_factor_a * parameters.replication_factor_b * parameters.replication_factor_c; // A ^ 3
    
    GuestParticle * guests = (GuestParticle *) malloc(MAX_GUESTS * sizeof(GuestParticle));
    int * N_g = (int *) calloc(2,sizeof(int)); // initialize current number of guests
    N_g[0] = 0; N_g[1] = 0;
    int N_g_total = 0; // total # guests
//  outputfile << "\tSize of guests = " << MAX_GUESTS * sizeof(particle_g) / (1024.0 * 1024.0) << " MB\n";
    
    //
    // Set up random number generators
    //
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    std::uniform_real_distribution<double> uniform01(0.0, 1.0); // uniformly distributed real no in [0,1]
    std::uniform_int_distribution<int> uniformint(0,10); // initialize with bogus range, will change later.
    std::uniform_int_distribution<int> type_generator(0, parameters.numadsorbates - 1); // for picking a component
    if (parameters.verbose) printf("Random number generator initialized...\n");

    //
    // Write settingsto outputfile
    //
    FILE * outputfile;
    char outputfilename[512];
    if (parameters.numadsorbates == 1)
        sprintf(outputfilename, "output_files/%s_%s_%sPa_gcmc.out", 
            parameters.frameworkname.c_str(), parameters.adsorbate[0].c_str(), argv[3]);
    if (parameters.numadsorbates == 2)
        sprintf(outputfilename, "output_files/%s_%s_%sPa_%s_%sPa_gcmc.out", 
            parameters.frameworkname.c_str(), 
            parameters.adsorbate[0].c_str(), argv[3],
            parameters.adsorbate[1].c_str(), argv[5]);
    outputfile = fopen(outputfilename, "w");
    WriteSettingsToOutputfile(outputfile, parameters, framework, forcefield, grid_info);
    if (parameters.verbose) printf("Wrote info to outputfile\n");

    //
    // Initialize adsorbate positions file to dump adsorbate positions
    //
    FILE * adsorbatepositionfile;
    char positionfilename[512];
    sprintf(positionfilename, "output_files/adsorbate_positions_%s.xyz", framework.name.c_str());
    adsorbatepositionfile = fopen(positionfilename, "w");
    int N_snapshots = 0;

    //
    //  Run GCMC
    //
    double start_of_sim_time=ReadTimer();
    if (parameters.verbose) printf("Starting simulation...\n");
    double E_gg_this_cycle = 0.0; double E_gf_this_cycle = 0.0; // assumes that we start with zero particles
    int cycle_counter = 0;
    int adsorbate_index_list[2][MAX_GUESTS] = { -1 }; // keep indices of particles here, initialze as -1
    
    for (int cycle = 0; cycle < parameters.numtrials; cycle++) {
        // each cycle corresponds to a number of Markov chain moves defined by ninnercycles
        int ninnercycles = (N_g[0] + N_g[1]) > 20 ? (N_g[0] + N_g[1]) : 20; // proportional to number of guests
        for (int inner_cycle = 0; inner_cycle < ninnercycles; inner_cycle ++) {
//          assert(N_g_total == N_g[0] + N_g[1]);
            cycle_counter += 1;

            double which_move = uniform01(generator); // particle swap or translation?
            double rand_for_acceptance = uniform01(generator); // for testing acceptance
            int which_type = type_generator(generator); // select adsorbate type 
            //
            //  MC trail: Insertion
            //
            if (which_move < parameters.p_exchange / 2.0) {
                stats.N_insertion_trials += 1;
                
                // insertion location @ these fractional coordinates
                guests[N_g_total].x_f = uniform01(generator) * parameters.replication_factor_a;
                guests[N_g_total].y_f = uniform01(generator) * parameters.replication_factor_b;
                guests[N_g_total].z_f = uniform01(generator) * parameters.replication_factor_c;
                // what are the Cartesian coords?
                double x, y, z; FractionalToCartesian(parameters.t_matrix, guests[N_g_total].x_f, guests[N_g_total].y_f, guests[N_g_total].z_f, x, y, z);
                guests[N_g_total].x = x;
                guests[N_g_total].y = y;
                guests[N_g_total].z = z;
                // declare particle type
                guests[N_g_total].type = which_type;

                // compute energy of this inserted particle
                double energy_gf = 1e6;
                if (which_type == 0)
                    energy_gf = GuestFrameworkEnergy(WrapToInterval(guests[N_g_total].x_f , 1.0),WrapToInterval(guests[N_g_total].y_f , 1.0), WrapToInterval(guests[N_g_total].z_f ,1.0),
                                    grid_info, energy_grid0);
                if (which_type == 1)
                    energy_gf = GuestFrameworkEnergy(WrapToInterval(guests[N_g_total].x_f , 1.0),WrapToInterval(guests[N_g_total].y_f , 1.0), WrapToInterval(guests[N_g_total].z_f ,1.0),
                                    grid_info, energy_grid1);
                double energy_gg = GuestGuestEnergy(N_g_total, N_g_total, guests, parameters);
                double E_insertion = energy_gf + energy_gg;

                // accept, loosely, if energetically favorable
                double acceptance_insertion = parameters.fugacity[which_type] * volume / ((N_g[which_type] + 1) * 1.3806488e7 * parameters.T) * exp(-E_insertion / parameters.T);
                if (parameters.debugmode) {
                    std::cout << "Adsorbate id list: ";
                    for (int ad = 0 ; ad < N_g[0] ; ad++)
                         std::cout << adsorbate_index_list[0][ad] << std::endl;
                    printf("\tProposal to insert at xf=%f,yf=%f,zf=%f\n.",guests[N_g_total].x_f,guests[N_g_total].y_f,guests[N_g_total].z_f);
                    std::cout << "\tE_gg = " << energy_gg << std::endl;
                    std::cout << "\tE_gf = " << energy_gf << std::endl;
                    std::cout << "\tAcceptance prob = " << acceptance_insertion << std::endl;
                }
                if (rand_for_acceptance < acceptance_insertion) {  // accept insertion move
                    if (parameters.debugmode) std::cout << "\tInsertion accepted. " << std::endl;
                    stats.N_insertions += 1;
                    if (energy_gf > 1e6)
                        std::cout << "Insertion accepted with huge energy" << std::endl;
                    // add adsorbate guests index to adsorbate_index_list
                    adsorbate_index_list[which_type][N_g[which_type]] = N_g_total;
                    N_g[which_type] += 1;
                    N_g_total += 1;
                    E_gg_this_cycle += energy_gg; E_gf_this_cycle += energy_gf;
                }
            }
            //
            //  MC trial: deletion 
            //
            else if (which_move < parameters.p_exchange)
            {
                if (parameters.debugmode) std::cout << "Deletion Trial." << std::endl;

                stats.N_deletion_trials += 1;
                if (N_g[which_type] > 0) {
                    // set new range for uniform int generator.
                    decltype(uniformint.param()) new_range(0, N_g[which_type] - 1);
                    uniformint.param(new_range);

                    // select idx of guest to delete in adsorbate_index_list
                    int idx_delete_type = uniformint(generator); 
//                  assert(idx_delete_type < N_g[which_type]);
                    int idx_delete = adsorbate_index_list[which_type][idx_delete_type]; // corresponding global ID in guests
                    if (idx_delete >= N_g_total) {
                        printf("Idx delete %d, N_g %d, Ng total %d,idx delete type %d\n",idx_delete,N_g[which_type],N_g_total,idx_delete_type);
                        printf("adsorbate index list: ");
                        for (int ad = 0 ; ad < N_g[0] ; ad++)
                         std::cout << adsorbate_index_list[0][ad] << std::endl;
                    }
//                  assert(idx_delete < N_g_total);

                    // compute energy of this guest that we propose to delete
                    double energy_gf = 1e6;
                    if (which_type == 0)
                        energy_gf = GuestFrameworkEnergy(WrapToInterval(guests[idx_delete].x_f , 1.0),WrapToInterval(guests[idx_delete].y_f , 1.0), WrapToInterval(guests[idx_delete].z_f ,1.0),
                                        grid_info, energy_grid0);
                    if (which_type == 1)
                        energy_gf = GuestFrameworkEnergy(WrapToInterval(guests[idx_delete].x_f , 1.0),WrapToInterval(guests[idx_delete].y_f , 1.0), WrapToInterval(guests[idx_delete].z_f ,1.0),
                                        grid_info, energy_grid1);
                    double energy_gg = GuestGuestEnergy(N_g_total, idx_delete, guests, parameters);
                    double E_deletion = energy_gf + energy_gg;

                    // accept deletion if, loosely, energetically favorable
                    double acceptance_del = (N_g[which_type] * 1.3806488e7 * parameters.T) / (parameters.fugacity[which_type] * volume) * exp(E_deletion / parameters.T);
                    if (rand_for_acceptance < acceptance_del) {
                        if (energy_gf > 1e6) {
                            std::cout << "Deletion accepted with huge energy" << std::endl;
                            printf("N_g = %d, energy_gg = %f, idx_delete = %d, idx_delete_type = %d, N_g1 = %d\n", N_g[0],energy_gg,idx_delete, idx_delete_type,N_g[1]);
                        }
                        stats.N_deletions += 1;
                        
                        // move last particle to here
                        guests[idx_delete] = guests[N_g_total - 1];
                        // if delete adsorbate in the middle of adsorbate_index_list, move the last one here so we see it.
                        adsorbate_index_list[which_type][idx_delete_type] = adsorbate_index_list[which_type][N_g[which_type] - 1];
                        // also we moved the last guest to the middle... find this!
                        for (int orz = 0 ; orz < N_g[guests[N_g_total - 1].type] ; orz ++) {
                            if (adsorbate_index_list[guests[N_g_total-1].type][orz] == N_g_total - 1) {
                                adsorbate_index_list[guests[N_g_total-1].type][orz] = idx_delete; //now this one is at idx_delete
                                break;
                            }
                        }

                        N_g[which_type] -= 1;
                        N_g_total -= 1;
                        E_gg_this_cycle -= energy_gg; E_gf_this_cycle -= energy_gf;
                        if (parameters.debugmode) {
                            std::cout << "Deletion accepted with probability " << acceptance_del << std::endl;
                        }
                    }
                }
            }
            //
            //  MC trial:  Translation
            //
            else if (which_move < parameters.p_exchange + parameters.p_move) {
                if (parameters.debugmode) 
                    printf("Translation Trial.\n");
                stats.N_move_trials += 1;

                if (N_g[which_type] > 0) {
                    // Randomly choose an adsorbate of which_type
                    decltype(uniformint.param()) new_range(0, N_g[which_type] - 1); // set new range for rng
                    uniformint.param(new_range);
                    int idx_move_type = uniformint(generator);
//                  assert(idx_move_type < N_g[which_type]);
                    int idx_move = adsorbate_index_list[which_type][idx_move_type]; // global ID in guests
//                  assert(idx_move < N_g_total);
                    
                    // compute energy in current (old) position
                    double energy_gf_old = 0.0;
                    if (which_type == 0)
                        energy_gf_old = GuestFrameworkEnergy(WrapToInterval(guests[idx_move].x_f , 1.0),WrapToInterval(guests[idx_move].y_f , 1.0), WrapToInterval(guests[idx_move].z_f ,1.0),
                                            grid_info, energy_grid0);
                    if (which_type == 1)
                        energy_gf_old = GuestFrameworkEnergy(WrapToInterval(guests[idx_move].x_f , 1.0),WrapToInterval(guests[idx_move].y_f , 1.0), WrapToInterval(guests[idx_move].z_f ,1.0),
                                            grid_info, energy_grid1);
                    double energy_gg_old = GuestGuestEnergy(N_g_total, idx_move, guests, parameters);

                    // store old coords
                    GuestParticle old_position = guests[idx_move];

                    // perturb the position in Cartesian space
                    guests[idx_move].x += parameters.delta * (uniform01(generator) - 0.5);
                    guests[idx_move].y += parameters.delta * (uniform01(generator) - 0.5);
                    guests[idx_move].z += parameters.delta * (uniform01(generator) - 0.5);

                    // convert back to fractional coords
                    double x_f, y_f, z_f;
                    CartesianToFractional(parameters.inv_t_matrix, x_f, y_f, z_f, guests[idx_move].x, guests[idx_move].y, guests[idx_move].z);
                    guests[idx_move].x_f = x_f;
                    guests[idx_move].y_f = y_f;
                    guests[idx_move].z_f = z_f;

                    // if move outside of box, reflect particle to the other side
                    if ((x_f > parameters.replication_factor_a) | (y_f > parameters.replication_factor_b) | (z_f > parameters.replication_factor_c) | 
                                        (x_f < 0.0) | (y_f < 0.0) | (z_f < 0.0)) {
                        // adjust fractional coordinates, then later Cartesian coords
                        guests[idx_move].x_f = WrapToInterval(x_f,parameters.replication_factor_a);
                        guests[idx_move].y_f = WrapToInterval(y_f,parameters.replication_factor_b);
                        guests[idx_move].z_f = WrapToInterval(z_f,parameters.replication_factor_c);
                        double x,y,z;
                        FractionalToCartesian(parameters.t_matrix, guests[idx_move].x_f, guests[idx_move].y_f, guests[idx_move].z_f, x, y, z);
                        guests[idx_move].x = x;
                        guests[idx_move].y = y;
                        guests[idx_move].z = z;
                    }
//                  assert(guests[idx_move].x_f > 0.0); assert(guests[idx_move].x_f < parameters.replication_factor_a);
//                  assert(guests[idx_move].y_f > 0.0); assert(guests[idx_move].y_f < parameters.replication_factor_b);
//                  assert(guests[idx_move].z_f > 0.0); assert(guests[idx_move].z_f < parameters.replication_factor_c);

                    // get energy at new position
                    double energy_gf_new = 1e6;
                    if (which_type == 0)
                        energy_gf_new = GuestFrameworkEnergy(WrapToInterval(guests[idx_move].x_f , 1.0),WrapToInterval(guests[idx_move].y_f , 1.0), WrapToInterval(guests[idx_move].z_f ,1.0),
                                            grid_info, energy_grid0);
                    if (which_type == 1)
                        energy_gf_new = GuestFrameworkEnergy(WrapToInterval(guests[idx_move].x_f , 1.0),WrapToInterval(guests[idx_move].y_f , 1.0), WrapToInterval(guests[idx_move].z_f ,1.0),
                                            grid_info, energy_grid1);
                    double energy_gg_new = GuestGuestEnergy(N_g_total, idx_move, guests, parameters);

                    double dE = energy_gf_new + energy_gg_new - (energy_gf_old + energy_gg_old);
                    if (rand_for_acceptance < exp(-dE/parameters.T)) {
                        stats.N_moves += 1; 
                        E_gg_this_cycle += energy_gg_new - energy_gg_old; E_gf_this_cycle += energy_gf_new - energy_gf_old;
                        if (energy_gf_new > 1e6) {
                            std::cout << "Move accepted with huge energy" << std::endl;
                            std::cout << "Old energy = " << energy_gf_old << std::endl;
                            std::cout << "New energy = " << energy_gf_new << std::endl;
                        }
                        // already overwrote coords with new coords, no need to update coords in this case
                    }
                    else {
                        // replace new cords with old coords
                        guests[idx_move] = old_position;
                    }
                } // end if N_g == 0
            } // end translation

            //
            //  PARTICLE IDENTITY CHANGE FOR DUAL COMPONENT
            //
            else if (which_move < parameters.p_exchange + parameters.p_move + parameters.p_identity_change) {
                if (N_g[which_type] > 0) { // if there are paricles of this type in the system
                    stats.N_ID_swap_trials += 1;

                    // pick which particles to change identity of
                    decltype(uniformint.param()) new_range (0, N_g[which_type] - 1); // set new range for rng
                    uniformint.param(new_range);
                    int idx_type = uniformint(generator); // which of this component?
                    int idx = adsorbate_index_list[which_type][idx_type]; // global index of this component
//                  assert(guests[idx].type == which_type);

                    // get energy difference with current ID
                    double energy_gg_old = GuestGuestEnergy(N_g_total, idx , guests, parameters);
                    double energy_gf_0 = GuestFrameworkEnergy(WrapToInterval(guests[idx].x_f , 1.0),WrapToInterval(guests[idx].y_f , 1.0), WrapToInterval(guests[idx].z_f ,1.0),
                                            grid_info, energy_grid0);
                    double energy_gf_1 = GuestFrameworkEnergy(WrapToInterval(guests[idx].x_f , 1.0),WrapToInterval(guests[idx].y_f , 1.0), WrapToInterval(guests[idx].z_f ,1.0),
                                            grid_info, energy_grid1);
                    
                    // change identity, recalculate energy
                    if (which_type == 0)
                        guests[idx].type = 1;
                    if (which_type == 1)
                        guests[idx].type = 0;
                    
                    double energy_gg_new = GuestGuestEnergy(N_g_total, idx , guests, parameters);
                    double dE;
                    if (which_type == 0) // if changing it to a 1
                        dE = energy_gf_1 + energy_gg_new - (energy_gf_0 + energy_gg_old);
                    if (which_type == 1) // if changing it to a 0
                        dE = energy_gf_0 + energy_gg_new - (energy_gf_1 + energy_gg_old);
                    int new_type = guests[idx].type;
                    double prob_acceptance_ID_swap = exp(-dE / parameters.T) * parameters.fugacity[new_type] / parameters.fugacity[which_type] * static_cast<double>(N_g[which_type]) / (static_cast<double>( N_g[new_type]) + 1.0);
                    if (rand_for_acceptance < prob_acceptance_ID_swap) { // ACCEPT
                        stats.N_ID_swaps += 1;
                        // keep guest particle as new declared type.
                        // keep track of energies
                        E_gg_this_cycle += energy_gg_new - energy_gg_old;
                        if (which_type == 0) // if changed to a 1
                            E_gf_this_cycle += energy_gf_1 - energy_gf_0;
                        if (which_type == 1) // if changed to a zero
                            E_gf_this_cycle += energy_gf_0 - energy_gf_1;
                        // deal with adorbate ID list
                        adsorbate_index_list[which_type][idx_type] = adsorbate_index_list[which_type][N_g[which_type] - 1]; // move last one of which_type to here, this one is deleted.
                        adsorbate_index_list[new_type][N_g[new_type]] = idx; // now this index is part of the other component ID
                        // update particle numbers
                        N_g[which_type] -= 1; // one less of which type
                        N_g[new_type] += 1; // one more of new type
                    }
                    else // if didn't accept
                        guests[idx].type = which_type; // go back to old type.
                } // end if N_g > 0
            } // end particle identity swap

            //
            //  MC trial: Regrow (take out molecule and insert it in a different, random location) like a move mostly
            //
            else {
                if (parameters.debugmode) 
                    printf("Regrow Trial.\n");
                stats.N_regrow_trials += 1;
                
                if (N_g[which_type] > 0) { // if there are paricles of this type in the system
                    // Randomly choose an adsorbate of which_type
                    decltype(uniformint.param()) new_range(0, N_g[which_type] - 1); // set new range for rng
                    uniformint.param(new_range);
                    int idx_move_type = uniformint(generator);
                    int idx_move = adsorbate_index_list[which_type][idx_move_type]; // global ID in guests particle array
                    
                    // compute energy in current (old) position
                    double energy_gf_old = 0.0;
                    if (which_type == 0)
                        energy_gf_old = GuestFrameworkEnergy(WrapToInterval(guests[idx_move].x_f , 1.0),
                                                             WrapToInterval(guests[idx_move].y_f , 1.0), 
                                                             WrapToInterval(guests[idx_move].z_f , 1.0),
                                                             grid_info, energy_grid0);
                    if (which_type == 1)
                        energy_gf_old = GuestFrameworkEnergy(WrapToInterval(guests[idx_move].x_f , 1.0),
                                                             WrapToInterval(guests[idx_move].y_f , 1.0), 
                                                             WrapToInterval(guests[idx_move].z_f , 1.0),
                                                             grid_info, energy_grid1);
                    double energy_gg_old = GuestGuestEnergy(N_g_total, idx_move, guests, parameters);
                    
                    // store old coords
                    GuestParticle old_position = guests[idx_move];

                    // insert this guest in a new random location
                    guests[N_g_total].x_f = uniform01(generator) * parameters.replication_factor_a;
                    guests[N_g_total].y_f = uniform01(generator) * parameters.replication_factor_b;
                    guests[N_g_total].z_f = uniform01(generator) * parameters.replication_factor_c;
                    // what are the Cartesian coords? assign them
                    double x, y, z; 
                    FractionalToCartesian(parameters.t_matrix, 
                                          guests[N_g_total].x_f, guests[N_g_total].y_f, guests[N_g_total].z_f, 
                                          x, y, z);
                    guests[N_g_total].x = x;
                    guests[N_g_total].y = y;
                    guests[N_g_total].z = z;
                    
                    // get energy at new position
                    double energy_gf_new = 1e6;
                    if (which_type == 0)
                        energy_gf_new = GuestFrameworkEnergy(WrapToInterval(guests[idx_move].x_f , 1.0),
                                                             WrapToInterval(guests[idx_move].y_f , 1.0), 
                                                             WrapToInterval(guests[idx_move].z_f , 1.0),
                                                             grid_info, energy_grid0);
                    if (which_type == 1)
                        energy_gf_new = GuestFrameworkEnergy(WrapToInterval(guests[idx_move].x_f , 1.0),
                                                             WrapToInterval(guests[idx_move].y_f , 1.0), 
                                                             WrapToInterval(guests[idx_move].z_f , 1.0),
                                                             grid_info, energy_grid1);
                    double energy_gg_new = GuestGuestEnergy(N_g_total, idx_move, guests, parameters);

                    // compute change in energy, accept move if "favorable"
                    double dE = energy_gf_new + energy_gg_new - (energy_gf_old + energy_gg_old);
                    if (rand_for_acceptance < exp(-dE/parameters.T)) {
                        stats.N_regrows += 1; 
                        E_gg_this_cycle += energy_gg_new - energy_gg_old; E_gf_this_cycle += energy_gf_new - energy_gf_old;
                        if (energy_gf_new > 1e6) {
                            std::cout << "Move accepted with huge energy" << std::endl;
                            std::cout << "Old energy = " << energy_gf_old << std::endl;
                            std::cout << "New energy = " << energy_gf_new << std::endl;
                        }
                        // already overwrote coords with new coords, no need to update coords in this case
                    }
                    else {
                        // replace new cords with old coords
                        guests[idx_move] = old_position;
                    }
                }
            }

            // assert N_g < MAX_GUESTS
            if (N_g_total >= MAX_GUESTS - 2) {
                printf("N_g > MAX_GUESTS!!!\n");
                exit(EXIT_FAILURE);
            }

            //
            // Collect statistics
            //
            if ((cycle > parameters.numequilibriumtrials) & (cycle_counter % parameters.samplefrequency == 0)) {
                stats.N_samples += 1;
                stats.N_g_avg[0] += N_g[0]; stats.N_g_avg[1] += N_g[1];
                stats.N_g2_avg[0] += N_g[0] * N_g[0]; stats.N_g2_avg[1] += N_g[1] * N_g[1];
                stats.guest_guest_energy_avg += E_gg_this_cycle; // divide by N_samples later
                stats.framework_guest_energy_avg += E_gf_this_cycle; // divide by N_samples later
            }
            
            //
            // Write adsorbate positions to file
            //
            if (parameters.writeadsorbatepositions) {
                if ((cycle > parameters.numequilibriumtrials) & (cycle_counter % parameters.writepositionfrequency == 0)) {
                    N_snapshots ++;
                    WriteGuestPostionsToFile(adsorbatepositionfile, guests, N_g_total, parameters);
                }
            }

            // for debug mode
            if (parameters.debugmode == true) {
                std::cout << "MC move " << cycle_counter << std::endl;
                std::cout << "\tN_g [0]= " << N_g[0] << std::endl;
                std::cout << "\tN_g [1]= " << N_g[1] << std::endl;
                std::cout << "\tN_g_total = " << N_g_total << std::endl;
                std::cout << "\tE_gg = " << E_gg_this_cycle << std::endl;
                std::cout << "\tE_gf = " << E_gf_this_cycle << std::endl;
            }


        } // end inner cycle

    }
    
    // take avg
    stats.guest_guest_energy_avg = 1.0 * stats.guest_guest_energy_avg / stats.N_samples;
    stats.framework_guest_energy_avg = 1.0 * stats.framework_guest_energy_avg / stats.N_samples;

    double sim_time = ReadTimer()-start_of_sim_time;
    fprintf(outputfile, "\nEnergy checks\n");
    // check energy calcs
    double E_gg_system = GuestGuestTotalEnergy(N_g[0] + N_g[1], guests, parameters);
    double E_gf_system = GuestFrameworkTotalEnergy(N_g[0] + N_g[1], energy_grid0, energy_grid1, guests, grid_info);
    fprintf(outputfile, "    E_gg total calc'ed at end: %f K\n", E_gg_system);
    fprintf(outputfile, "    E_gg from adding dE's throughout simulation: %f K\n", E_gg_this_cycle);
    fprintf(outputfile, "    E_gf total calc'ed at end: %f K\n", E_gf_system);
    fprintf(outputfile, "    E_gf from adding dE's throughout simulation: %f K\n", E_gf_this_cycle);
    
    // write stats
    fprintf(outputfile, "\nGCMC Statistics\n");
    fprintf(outputfile, "    Simulation time: %f min\n", sim_time / 60.0);
    fprintf(outputfile, "    Insertions: %d / %d\n", stats.N_insertions, stats.N_insertion_trials);
    fprintf(outputfile, "    Deletions: %d / %d\n", stats.N_deletions, stats.N_deletion_trials);
    fprintf(outputfile, "    Moves: %d / %d\n", stats.N_moves, stats.N_move_trials);
    fprintf(outputfile, "    ID swaps: %d / %d\n", stats.N_ID_swaps, stats.N_ID_swap_trials);
    fprintf(outputfile, "    Regrows: %d / %d\n", stats.N_regrows, stats.N_regrow_trials);
    fprintf(outputfile, "\n    Number of samples: %d\n", stats.N_samples);
    
    // write loadings
    for (int n_c = 0; n_c < parameters.numadsorbates; n_c ++) {
        stats.N_g_avg[n_c] = 1.0 * stats.N_g_avg[n_c] / stats.N_samples;
        stats.N_g2_avg[n_c] = 1.0 * stats.N_g2_avg[n_c] / stats.N_samples;
        double N_confidence_bound = sqrt(stats.N_g2_avg[n_c] - stats.N_g_avg[n_c] * stats.N_g_avg[n_c])/sqrt(1.0 * stats.N_samples); // sigma / sqrt(N)

        fprintf(outputfile, "\n    Adsorbate: %s\n", parameters.adsorbate[n_c].c_str());
        fprintf(outputfile, "        Fugacity = %f Pa\n", parameters.fugacity[n_c]);
        fprintf(outputfile, "        <N_g> (%s) = %f +/- %f molecules/ unit cell\n", parameters.adsorbate[n_c].c_str(), 1.0 * stats.N_g_avg[n_c] / parameters.replication_factor_a / parameters.replication_factor_b / parameters.replication_factor_c, N_confidence_bound);
        fprintf(outputfile, "        <N_g> (%s) = %f moles/m3\n", parameters.adsorbate[n_c].c_str(), stats.N_g_avg[n_c] / volume / 6.022e-7);
        fprintf(outputfile, "        <N_g> (%s) = %f moles/kg framework\n", parameters.adsorbate[n_c].c_str(), stats.N_g_avg[n_c]/volume/6.022e-7/framework.density);
    }
    fprintf(outputfile, "\n     <E_gg> = %f kJ/mol = %f K\n", stats.guest_guest_energy_avg * 8.314 / 1000.0, stats.guest_guest_energy_avg);
    fprintf(outputfile, "     <E_gf> = %f kJ/mol = %f K", stats.framework_guest_energy_avg * 8.314 / 1000.0, stats.framework_guest_energy_avg);

    if (parameters.writeadsorbatepositions) 
        fprintf(outputfile, "\nWrote adsorbate positions in xyz file every %d MC moves\n    Took total of %d snapshots\n", parameters.writepositionfrequency, N_snapshots);

    fclose(outputfile); 
}
