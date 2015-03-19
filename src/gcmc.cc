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
// TODO difference between erbose and debug?

void InitializeGCMCStats(GCMCStats & stats) {
    // initializes GCMC statistics to zero
    stats.N_move_trials = 0; stats.N_insertion_trials = 0; stats.N_deletion_trials = 0; stats.N_ID_swap_trials = 0; stats.N_regrow_trials = 0;
    stats.N_insertions = 0; stats.N_deletions = 0; stats.N_moves = 0; stats.N_ID_swaps = 0; stats.N_regrows = 0;
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

bool OutsideUnitCell(double x_f, double y_f, double z_f,
                    GCMCParameters parameters) {
    // check if outside bounding box
    if ((x_f > parameters.replication_factor_a) | 
        (y_f > parameters.replication_factor_b) | 
        (z_f > parameters.replication_factor_c) | 
        (x_f < 0.0) | (y_f < 0.0) | (z_f < 0.0)) {
        return true;
    }
    else
        return false;
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

//void DeleteGuest(GuestMolecule * guestmolecules,
//                 GuestBead * guestbeadsint,
//                 GuestMoleculeInfo * guestmoleculeinfo,
//                 int guestmoleculeid, 
//                 int N_g) {
//    // Remove guest guestmoleculeid
//    
//    // Remove each bead from the guest
//    for (int b = 0; b < guestmoleculeinfo[guestmolecules[guestmoleculeid]].nbeads; b++ ) {  // for each bead that makes up this guest molecule
//        int beadid = guestmolecules[guestmoleculeid].beadID[b];
//
//    }
//}

double GuestGuestEnergy(int N_g,
                        int guestmoleculeid, 
                        GuestMoleculeInfo * guestmoleculeinfo,
                        GuestMolecule * guestmolecules,
                        GuestBead * guestbeads,
                        GCMCParameters parameters) {
    // Compute potential energy of guest molecule guestmoleculeid, the contribution from other guests.
    double E_gg = 0.0; // initiate
    
    // for each bead in this guest molecule
    for (int b_this = 0; b_this < guestmoleculeinfo[guestmolecules[guestmoleculeid].type].nbeads; b_this++) {
        // coordinates of this bead
        double x_f = guestbeads[guestmolecules[guestmoleculeid].beadID[b_this]].x_f;
        double y_f = guestbeads[guestmolecules[guestmoleculeid].beadID[b_this]].y_f;
        double z_f = guestbeads[guestmolecules[guestmoleculeid].beadID[b_this]].z_f;
        // get bead type for this
        int this_bead_type = guestbeads[guestmolecules[guestmoleculeid].beadID[b_this]].type;
        // for each other guest molecule 
        for (int k = 0 ; k < N_g; k++) {
            // do not include self interation orz
            if (k == guestmoleculeid) 
                continue; 
            // get energy contribution from each bead in other guest molecule
            for (int b = 0; b < guestmoleculeinfo[guestmolecules[k].type].nbeads; b++ ) {
                // get other bead type for other guest molecule
                int other_bead_type = guestbeads[guestmolecules[k].beadID[b]].type;

                // distance in fractional coords
                double dx_f = guestbeads[guestmolecules[k].beadID[b]].x_f - x_f;
                double dy_f = guestbeads[guestmolecules[k].beadID[b]].y_f - y_f;
                double dz_f = guestbeads[guestmolecules[k].beadID[b]].z_f - z_f;

                // take nearest image
                dx_f = dx_f - parameters.replication_factor_a * round(dx_f / parameters.replication_factor_a);
                dy_f = dy_f - parameters.replication_factor_b * round(dy_f / parameters.replication_factor_b);
                dz_f = dz_f - parameters.replication_factor_c * round(dz_f / parameters.replication_factor_c);
                // assert within bounds #todo remove later
                assert (dx_f < 0.5 * parameters.replication_factor_a);
                assert (dy_f < 0.5 * parameters.replication_factor_b);
                assert (dz_f < 0.5 * parameters.replication_factor_c);
                assert (dx_f > -0.5 * parameters.replication_factor_a);
                assert (dy_f > -0.5 * parameters.replication_factor_b);
                assert (dz_f > -0.5 * parameters.replication_factor_c);

                // distance in Cartesian
                double dx, dy, dz;
                FractionalToCartesian(parameters.t_matrix,
                                        dx_f, dy_f, dz_f,
                                        dx, dy, dz);
                
                // Compute LJ potential
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < min_r) 
                    return 100000000000000.0; // overwrite if too small
                if (r2 < parameters.r_cutoff_squared) {
                    double sigma_over_r_sixth = pow(parameters.sigma_squared_matrix[this_bead_type][other_bead_type] / r2, 3.0); 
                    E_gg += 4.0 * parameters.epsilon_matrix[this_bead_type][other_bead_type] * sigma_over_r_sixth * (sigma_over_r_sixth - 1.0);
                }
            }  // end loop over BEADS of other guest molecule
        }  // end loop over other guest molecules
    }  // end loop over beads of THIS guest molecule
    return E_gg;
}

double TotalGuestGuest(int N_g,
                        GuestMoleculeInfo * guestmoleculeinfo,
                        GuestMolecule * guestmolecules,
                        GuestBead * guestbeads,
                        GCMCParameters parameters) {
    // Calculate total system guest-guest energy
    double E_gg = 0.0;
    for (int i = 0; i < N_g; i ++) {
        E_gg += GuestGuestEnergy(N_g,
                        i, 
                        guestmoleculeinfo,
                        guestmolecules,
                        guestbeads,
                        parameters);
    }
    return E_gg / 2.0; // 2 is for double counting
}

int GetEnergyGridIndex(int i, int j, int k, GridInfo grid_info) {
    //
    //  Get index of flattened, 1D pointer array that corresponds to grid point index (i, j, k)
    //
    // returns index in energy_grid corresponding to gridpoint index i,j,k
    return k + j * grid_info.N_z + i * grid_info.N_y * grid_info.N_z;
}

//
// Interpolate energy grid to get bead-framework energy
//  
double BeadFrameworkEnergy(double x_f_, double y_f_, double z_f_,
            GridInfo grid_info,
            double * energy_grid) {
    // reflect fractional coords in [0,1] for grid
    // TODO: put wrap function here?
    assert ((x_f_ > 0.0) & (y_f_ > 0.0) & (z_f_ > 0.0));
    assert ((x_f_ < 1.0) & (y_f_ < 1.0) & (z_f_ < 1.0));
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

double GuestFrameworkEnergy(int guestmoleculeid, 
                            GuestMoleculeInfo * guestmoleculeinfo,
                            GuestMolecule * guestmolecules,
                            GuestBead * guestbeads,
                            GridInfo grid_info,
                            double ** energy_grids) {
    // Compute energy of guestmoleculeid'th guest molecule
    double E_gf = 0.0;
    for (int b = 0; b < guestmoleculeinfo[guestmolecules[guestmoleculeid].type].nbeads; b++) {  // for each bead that makes up this guest molecule
        int beadid = guestmolecules[guestmoleculeid].beadID[b];  // id of bead in guestbeads
        int beadtype = guestbeads[beadid].type;  // int id of bead type
        E_gf += BeadFrameworkEnergy(WrapToInterval(guestbeads[beadid].x_f, 1.0),
                                    WrapToInterval(guestbeads[beadid].y_f, 1.0),
                                    WrapToInterval(guestbeads[beadid].z_f, 1.0),
                                    grid_info, energy_grids[beadtype]);
    }
    return E_gf;
}

double TotalGuestFrameworkEnergy(int N_g, 
                            GuestMoleculeInfo * guestmoleculeinfo,
                            GuestMolecule * guestmolecules,
                            GuestBead * guestbeads,
                            GridInfo grid_info,
                            double ** energy_grids) {
    // Calculate total system framework-guest energy
    double E_gf = 0.0;
    for (int i = 0; i < N_g; i ++) {
        E_gf += GuestFrameworkEnergy(i,
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            grid_info,
                            energy_grids);
    }
    return E_gf;
}

//
// Write guest positions to file
//
void WriteGuestPostionsToFile(FILE * positionfile, 
                              int N_g_total,
                              GuestMoleculeInfo * guestmoleculeinfo,
                              GuestMolecule * guestmolecules,
                              GuestBead * guestbeads,
                              GCMCParameters parameters) {
    // for each guest molecule...
    for (int i = 0; i < N_g_total ; i++) {
        // for each bead in this guest molecule...
        for (int b = 0; b < guestmoleculeinfo[guestmolecules[i].type].nbeads; b++) {
            int beadid = guestmolecules[i].beadID[b];  // ID of this bead
            int beadtype = guestbeads[beadid].type;
            std::string beadlabel = (guestmoleculeinfo[guestmolecules[i].type].beadtypes[0] == beadtype) ? 
                                    guestmoleculeinfo[guestmolecules[i].type].beadlabels[0] : 
                                    guestmoleculeinfo[guestmolecules[i].type].beadlabels[1]; 
            fprintf(positionfile, "%s %f %f %f\n", 
                    beadlabel.c_str(),
                    guestbeads[beadid].x, guestbeads[beadid].y, guestbeads[beadid].z);
        }
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
    // Read adsorbate info
    //
    GuestMoleculeInfo guestmoleculeinfo[2];  // at most 2 adsorbates
    GetGuestMoleculeInfo(guestmoleculeinfo, parameters);
    if (parameters.verbose) printf("Loaded guest molecule information\n");
    
    //
    // Get adsorbate LJ params and store in epsilon and sigma_squared matrices
    //
    for (int bx = 0; bx < parameters.nuniquebeads; bx++) {
        PairDouble eps_sig_x = GrabGuestForceFieldParams(forcefield, parameters.uniquebeadlist[bx]); 
        for (int by = 0; by < parameters.nuniquebeads; by++) {
            PairDouble eps_sig_y = GrabGuestForceFieldParams(forcefield, parameters.uniquebeadlist[by]);
            parameters.epsilon_matrix[bx][by] = sqrt(eps_sig_x.arg1 * eps_sig_y.arg1); 
            double sigma_mixed = (eps_sig_x.arg2 + eps_sig_y.arg2) / 2.0;
            parameters.sigma_squared_matrix[bx][by] = sigma_mixed * sigma_mixed;
        }
    }
    if (parameters.verbose) printf("Fetched adsorbate FF parameters\n");

    //
    // Import energy grid
    // + pocket blocking if enabled.
    //
    GridInfo grid_info; // for storing grid info
    double ** energy_grids = (double **) malloc(parameters.nuniquebeads * sizeof(double *));
    bool pocket_block_verbose_mode = false; // writes grid before and after
    for (int n_c = 0; n_c < parameters.nuniquebeads; n_c ++) {
        int N_x_temp, N_y_temp, N_z_temp;
        if (parameters.verbose) printf("Importing energy grid %d\n", n_c);
        char gridfilename[512];
        sprintf(gridfilename, "data/grids/%s_%s_%s.txt", framework.name.c_str(), parameters.uniquebeadlist[n_c].c_str(), forcefield.name.c_str());

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
        else {
            // assert other energy grids hv same resolution...
            assert(N_x_temp = grid_info.N_x);
            assert(N_y_temp = grid_info.N_y);
            assert(N_z_temp = grid_info.N_z);
        }

        // create grid import format
        energy_grids[n_c] = (double *) calloc(grid_info.numtotalpts, sizeof(double));

        int i = 0; int j = 0;
        while(getline(gridfile,line)) {
            std::istringstream this_line(line);
            for (int k = 0; k < grid_info.N_z; k ++) {  // each line is a pencil of z's
                int index_here = k + j * grid_info.N_z + i * grid_info.N_y * grid_info.N_z;
                this_line >> energy_grids[n_c][index_here];
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
        if (parameters.verbose) 
            printf("energy grid %d for bead %s imported successfully.\n", n_c, parameters.uniquebeadlist[n_c].c_str());
        
        //
        // Flood fill/ pocket blocking, if enabled
        //
        if (parameters.pocketblocking) {
            double time_before = ReadTimer();
            if (pocket_block_verbose_mode) {
                printf("Pocket blocking beginning. Write a cube for before\n");
                if (n_c == 0)
                    WriteCube("before_blocking_0", framework, parameters, energy_grids[0], grid_info); //just so you can see the grid before ...
            }

            // pass energy grid to FindAndBlockPockets, which will do the job.
            grid_info.numpockets[n_c] = FindAndBlockPockets(energy_grids[n_c], grid_info, parameters.T, parameters);
            
            if (pocket_block_verbose_mode) {
                printf("Pocket blocking finished. Write a cube for after\n");
                double time_after = ReadTimer();
                printf("Time spent to find and block pockets: %f s\n", time_after - time_before);
                if (n_c == 0)
                    WriteCube("after_blocking_0", framework, parameters, energy_grids[0], grid_info); // ... and after pocket blocking
            }
        }
    }
    
    //
    // Initialize stats and guests array (includes both components via "type" attribute)
    //
    GCMCStats stats;
    InitializeGCMCStats(stats);
    GuestMolecule * guestmolecules = (GuestMolecule *) malloc(MAX_GUESTS * sizeof(GuestMolecule));
    GuestBead * guestbeads = (GuestBead *) malloc(MAX_GUESTS * 2 * sizeof(GuestBead));  // each guest can hv up to two beads
    int * N_g = (int *) calloc(2, sizeof(int)); // initialize current number of guests
    N_g[0] = 0; N_g[1] = 0;
    int N_beads = 0;  // number of beads
    int N_g_total = 0; // total # guests
    double volume = framework.volume_unitcell * parameters.replication_factor_a * parameters.replication_factor_b * parameters.replication_factor_c; // A ^ 3
//  outputfile << "\tSize of guests = " << MAX_GUESTS * sizeof(particle_g) / (1024.0 * 1024.0) << " MB\n";
    
    //
    // Set up random number generators
    //
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    std::uniform_real_distribution<double> uniform01(0.0, 1.0); // uniformly distributed real no in [0,1]
    std::uniform_int_distribution<int> uniformint(0, 10); // initialize with bogus range, will change later.
    std::uniform_int_distribution<int> beadpicker(0, 1); // For picking a random bead in particle ID exchange
    std::uniform_int_distribution<int> type_generator(0, parameters.numadsorbates - 1); // for picking a component
    if (parameters.verbose) 
        printf("Random number generator initialized...\n");

    //
    // Write settings to outputfile
    //
    FILE * outputfile;
    char outputfilename[512];
    if (parameters.numadsorbates == 1)
        sprintf(outputfilename, "output_files/%s_%s_%sPa_%gK_gcmc.out", 
            parameters.frameworkname.c_str(), parameters.adsorbate[0].c_str(), argv[3], parameters.T);
    if (parameters.numadsorbates == 2)
        sprintf(outputfilename, "output_files/%s_%s_%sPa_%s_%sPa_%gK_gcmc.out", 
            parameters.frameworkname.c_str(), 
            parameters.adsorbate[0].c_str(), argv[3],
            parameters.adsorbate[1].c_str(), argv[5],
            parameters.T);
    outputfile = fopen(outputfilename, "w");
    WriteSettingsToOutputfile(outputfile, parameters, framework, forcefield, grid_info, guestmoleculeinfo);
    if (parameters.verbose) printf("Wrote info to outputfile\n");

    //
    // Initialize adsorbate positions file to dump adsorbate positions
    //
    FILE * adsorbatepositionfile;
    int N_snapshots = 0;
    if (parameters.writeadsorbatepositions) {
        char positionfilename[1024];
        sprintf(positionfilename, "output_files/adsorbate_positions_%s.xyz", framework.name.c_str());
        adsorbatepositionfile = fopen(positionfilename, "w");
    }

    //
    //  Run GCMC
    //
    double start_of_sim_time = ReadTimer();
    if (parameters.verbose) 
        printf("Starting simulation...\n");
    double E_gg_this_cycle = 0.0; 
    double E_gf_this_cycle = 0.0; // assumes that we start with zero particles
    int cycle_counter = 0;

    // In this list, store index in guestmolecules that each guest type has
    int guestmolecule_index_list[2][MAX_GUESTS] = { -1 }; // keep indices of particles here, initialze as -1

//    printf("Guest framework energy (.2,.2,.3)=%f\n", 
//        BeadFrameworkEnergy(WrapToInterval(.2, 1.0),
//                                    WrapToInterval(.2, 1.0),
//                                    WrapToInterval(.2, 1.0),
//                                    grid_info, energy_grids[1]));
    
    for (int cycle = 0; cycle < parameters.numtrials; cycle++) {
        // each cycle corresponds to a number of Markov chain moves defined by ninnercycles
        int ninnercycles = N_g_total > 20 ? N_g_total : 20; // proportional to number of guests
        for (int inner_cycle = 0; inner_cycle < ninnercycles; inner_cycle ++) {
//          assert(N_g_total == N_g[0] + N_g[1]);
            cycle_counter += 1;
            
            // generate random numbers
            double whichMCmove = uniform01(generator); // Insertion, deletion, translation, ID swap, or regrow?
            double rand_for_acceptance = uniform01(generator); // for testing acceptance
            int which_type = type_generator(generator); // select guest type 
            double theta, phi;  // angles on sphere for two bead molecule
            if (guestmoleculeinfo[which_type].nbeads > 1) {
                // generate randomly distributed point on a sphere
                // http://mathworld.wolfram.com/SpherePointPicking.html
                theta = 2 * M_PI * uniform01(generator);
                phi = acos(2 * uniform01(generator) - 1);
            }
            
            //
            //  MC trial: Insertion
            //
            if (whichMCmove < parameters.p_exchange / 2.0) {
                stats.N_insertion_trials += 1;

                // add new guest, declare particle type
                guestmolecules[N_g_total].type = which_type;
                
                // add this bead to guests
                guestmolecules[N_g_total].beadID[0] = N_beads;
                // insert first bead @ these fractional coordinates
                guestbeads[N_beads].x_f = uniform01(generator) * parameters.replication_factor_a;
                guestbeads[N_beads].y_f = uniform01(generator) * parameters.replication_factor_b;
                guestbeads[N_beads].z_f = uniform01(generator) * parameters.replication_factor_c;
                // what are the Cartesian coords? updat guests array
                double x, y, z; 
                FractionalToCartesian(parameters.t_matrix, 
                                      guestbeads[N_beads].x_f, guestbeads[N_beads].y_f, guestbeads[N_beads].z_f, 
                                      x, y, z);
                guestbeads[N_beads].x = x; 
                guestbeads[N_beads].y = y; 
                guestbeads[N_beads].z = z; 
                // define which guest this bead belongs to
                guestbeads[N_beads].guestmoleculeID = N_g_total;
                // define bead type for energy computations
                guestbeads[N_beads].type = guestmoleculeinfo[which_type].beadtypes[0];

                // add second bead if not LJ sphere
                if (guestmoleculeinfo[which_type].nbeads > 1) {
                    // generate Cartesian coords of second bead
                    guestbeads[N_beads + 1].x = guestbeads[N_beads].x + guestmoleculeinfo[which_type].bondlength * sin(phi) * cos(theta);
                    guestbeads[N_beads + 1].y = guestbeads[N_beads].y + guestmoleculeinfo[which_type].bondlength * sin(phi) * sin(theta);
                    guestbeads[N_beads + 1].z = guestbeads[N_beads].z + guestmoleculeinfo[which_type].bondlength * cos(phi);

                    // convert to fractional 
                    double x_f, y_f, z_f;
                    CartesianToFractional(parameters.inv_t_matrix, 
                                          x_f, y_f, z_f, 
                                          guestbeads[N_beads + 1].x, guestbeads[N_beads + 1].y, guestbeads[N_beads + 1].z);
                    guestbeads[N_beads + 1].x_f = x_f;
                    guestbeads[N_beads + 1].y_f = y_f;
                    guestbeads[N_beads + 1].z_f = z_f;  // TODO pass this directly to cartesiantofrac
                    
                    // SECOND BEAD IS ALLOWED OUTSIDE OF THE UNIT CELL! (not the first though) 
                     
                    // define which guest this bead belongs to
                    guestbeads[N_beads + 1].guestmoleculeID = N_g_total;
                    // define bead type for energy computations
                    guestbeads[N_beads + 1].type = guestmoleculeinfo[which_type].beadtypes[1];
                    
                    // add this bead to guests
                    guestmolecules[N_g_total].beadID[1] = N_beads + 1;
                }

                // compute energy of this inserted particle
                double E_gf = GuestFrameworkEnergy(N_g_total, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            grid_info,
                            energy_grids);  // framework-guest
                double E_gg = GuestGuestEnergy(N_g_total,
                        N_g_total, 
                        guestmoleculeinfo,
                        guestmolecules,
                        guestbeads,
                        parameters);  // guest-guest
                double E_insertion = E_gf + E_gg;
                
                // for debug mode, print stuff off
                if (parameters.debugmode) {
                    printf("INSERTION PROPOSAL of type %d\n", which_type);
                    printf("\n\tFirst bead insertion proposal at xf=%f,yf=%f,zf=%f\n", 
                                guestbeads[N_beads].x_f, guestbeads[N_beads].y_f, guestbeads[N_beads].z_f);
                    if (guestmoleculeinfo[which_type].nbeads > 1) {
                        printf("\tSecond bead at xf=%f,yf=%f,zf=%f\n.", 
                                    guestbeads[N_beads+1].x_f, guestbeads[N_beads+1].y_f, guestbeads[N_beads+1].z_f);
                        printf("\tBond length = %f\n", sqrt( pow(guestbeads[N_beads].x - guestbeads[N_beads+1].x, 2) +
                                                             pow(guestbeads[N_beads].y - guestbeads[N_beads+1].y, 2) +
                                                             pow(guestbeads[N_beads].z - guestbeads[N_beads+1].z, 2)));

                    }
                    std::cout << "\tE_gg = " << E_gg << std::endl;
                    std::cout << "\tE_gf = " << E_gf << std::endl;
                }
                
                // accept, loosely, if energetically favorable
                double acceptance_insertion = parameters.fugacity[which_type] * volume / ((N_g[which_type] + 1) * 1.3806488e7 * parameters.T) * exp(-E_insertion / parameters.T);
                if (parameters.debugmode) 
                    std::cout << "\tAcceptance prob = " << acceptance_insertion << std::endl;
                if (rand_for_acceptance < acceptance_insertion) {  
                    if (parameters.debugmode) 
                        printf("\tInsertion ACCEPTED.\n");

                    stats.N_insertions += 1;
                    if (E_insertion > 1e6)
                        std::cout << "Insertion accepted with huge energy" << std::endl;
                    // add adsorbate guests index to guestmolecule_index_list
                    guestmolecule_index_list[which_type][N_g[which_type]] = N_g_total;
                    // update molecule and bead count
                    N_g[which_type] += 1;
                    N_g_total += 1;
                    N_beads += guestmoleculeinfo[which_type].nbeads;
                    // update system energies
                    E_gg_this_cycle += E_gg; 
                    E_gf_this_cycle += E_gf;
                }
            }

            //
            //  MC trial: deletion 
            //
            else if (whichMCmove < parameters.p_exchange)
            {
                if (parameters.debugmode) 
                    printf("DELETION PROPOSAL of type %d\n", which_type);
                stats.N_deletion_trials += 1;
                
                if (N_g[which_type] > 0) {
                    // set new range for uniform int generator [0 - # of this type of guest - 1]
                    decltype(uniformint.param()) new_range(0, N_g[which_type] - 1);
                    uniformint.param(new_range);

                    // randomly select guest of this type to propose to delete
                    int idx_thisguesttype = uniformint(generator);  // index in guestmolecule_index_list
                    int idx_guestmolecule = guestmolecule_index_list[which_type][idx_thisguesttype]; // corresponding global ID in guests
                    assert(idx_guestmolecule < N_g_total);
                    if (parameters.debugmode)
                        printf("\tProposal to delete guest molecule %d\n", idx_guestmolecule);
                    
                    // compute energy of this guest that we propose to delete
                    double E_gf = GuestFrameworkEnergy(idx_guestmolecule, 
                                guestmoleculeinfo,
                                guestmolecules,
                                guestbeads,
                                grid_info,
                                energy_grids);  // framework-guest
                    double E_gg = GuestGuestEnergy(N_g_total,
                            idx_guestmolecule, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            parameters);  // guest-guest
                    double E_deletion = E_gf + E_gg;
                    
                    // accept deletion if, loosely, energetically favorable
                    double acceptance_del = (N_g[which_type] * 1.3806488e7 * parameters.T) / (parameters.fugacity[which_type] * volume) * exp(E_deletion / parameters.T);
                    if (rand_for_acceptance < acceptance_del) {
                        stats.N_deletions += 1;
                        
                        if (E_deletion > 1e6) {
                            std::cout << "Deletion accepted with huge energy" << std::endl;
                            printf("N_g = %d, energy_gg = %f, idx_delete = %d, idx_delete_type = %d, N_g1 = %d\n", N_g[0], E_gg,idx_guestmolecule, idx_thisguesttype, N_g[1]);
                        }
                        
                        // replace this deleted guest's beads with the beads at the end of guestbeads
                        for (int b = 0; b < guestmoleculeinfo[which_type].nbeads; b++) {
                            int beadid = guestmolecules[idx_guestmolecule].beadID[b];
                            // move last bead to here
                            guestbeads[beadid] = guestbeads[N_beads - 1 - b];
                            // change molecule ID of beads for this guest
                            int otherguestmoleculeid =  guestbeads[N_beads - 1 - b].guestmoleculeID;
                            // change beadID in this other guest molecule to beadID
                            bool found = false;
                            for (int bi = 0; bi < guestmoleculeinfo[guestmolecules[otherguestmoleculeid].type].nbeads; bi++) {
                                if (guestmolecules[otherguestmoleculeid].beadID[bi] == N_beads - 1 - b) {
                                    guestmolecules[otherguestmoleculeid].beadID[bi] = beadid;
                                    found = true;
                                }
                            }
                            assert(found == true); // TODO remove later, for checking
                        }
                        // move last guest molecule to here
                        guestmolecules[idx_guestmolecule] = guestmolecules[N_g_total - 1];
                        // change guestmoleculeID for its beads
                        for (int b = 0; b < guestmoleculeinfo[guestmolecules[N_g_total - 1].type].nbeads; b++) 
                                guestbeads[guestmolecules[N_g_total - 1].beadID[b]].guestmoleculeID = idx_guestmolecule;
                        
                        // if delete adsorbate in the middle of guestmolecule_index_list, move the last one here so we see it.
                        guestmolecule_index_list[which_type][idx_thisguesttype] = guestmolecule_index_list[which_type][N_g[which_type] - 1];
                        // also we moved the last guest to the middle... find this!
                        for (int orz = 0 ; orz < N_g[guestmolecules[N_g_total - 1].type]; orz ++) {
                            if (guestmolecule_index_list[guestmolecules[N_g_total-1].type][orz] == N_g_total - 1) {
                                //now this one is at idx_guestmolecule
                                guestmolecule_index_list[guestmolecules[N_g_total-1].type][orz] = idx_guestmolecule; 
                                break;
                            }
                        }
                        
                        // update guest and bead counts
                        N_g[which_type] -= 1;
                        N_g_total -= 1;
                        N_beads -= guestmoleculeinfo[guestmolecules[idx_guestmolecule].type].nbeads;
                        // update system energies
                        E_gg_this_cycle -= E_gg; 
                        E_gf_this_cycle -= E_gf;
                        if (parameters.debugmode) {
                            std::cout << "Deletion accepted with probability " << acceptance_del << std::endl;
                        }
                    }
                }
            }
            
            //
            //  MC trial:  Translation
            //
            else if (whichMCmove < parameters.p_exchange + parameters.p_move) {
                if (parameters.debugmode) 
                    printf("Translation Trial.\n");
                stats.N_move_trials += 1;

                if (N_g[which_type] > 0) {
                    // Randomly choose an adsorbate of which_type
                    decltype(uniformint.param()) new_range(0, N_g[which_type] - 1); // set new range for rng
                    uniformint.param(new_range);
                    int idx_move_type = uniformint(generator);  // corresponds to index of this type in guestmolecule_index
                    int idx_guestmolecule = guestmolecule_index_list[which_type][idx_move_type]; // global ID in guests
                    
                    // compute energy in current (old) position
                    double E_gf_old = GuestFrameworkEnergy(idx_guestmolecule, 
                                guestmoleculeinfo,
                                guestmolecules,
                                guestbeads,
                                grid_info,
                                energy_grids);  // framework-guest
                    double E_gg_old = GuestGuestEnergy(N_g_total,
                            idx_guestmolecule, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            parameters);  // guest-guest
                    double E_old = E_gf_old + E_gg_old;

                    // get perturbation of the position in Cartesian space
                    double x_perturb = parameters.delta * (uniform01(generator) - 0.5);
                    double y_perturb = parameters.delta * (uniform01(generator) - 0.5);
                    double z_perturb = parameters.delta * (uniform01(generator) - 0.5);

                    GuestBead old_beads[2];  // for old coords of beads
                    // move each bead of this guest the same amount
                    for (int b = 0; b < guestmoleculeinfo[which_type].nbeads; b++) {
                        int beadid = guestmolecules[idx_guestmolecule].beadID[b];
                        // store old coords of beads
                        old_beads[b] = guestbeads[beadid];
                        
                        // perturb bead
                        guestbeads[beadid].x += x_perturb;
                        guestbeads[beadid].y += y_perturb;
                        guestbeads[beadid].z += z_perturb;
                        // convert back to fractional coords
                        double x_f, y_f, z_f;
                        CartesianToFractional(parameters.inv_t_matrix, 
                                             x_f, y_f, z_f, 
                                             guestbeads[beadid].x, guestbeads[beadid].y, guestbeads[beadid].z);
                        guestbeads[beadid].x_f = x_f;
                        guestbeads[beadid].y_f = y_f;
                        guestbeads[beadid].z_f = z_f;

                    }
                    
                    // if FIRST BEAD is moved outside of box, reflect ENTIRE guest to the other side
                    int first_beadid = guestmolecules[idx_guestmolecule].beadID[0];
                    double x_f_correction = 0.0; double y_f_correction = 0.0; double z_f_correction = 0.0;
                    if (guestbeads[first_beadid].x_f > 1.0 * parameters.replication_factor_a)
                        x_f_correction = - 1.0 * parameters.replication_factor_a; 
                    if (guestbeads[first_beadid].y_f > 1.0 * parameters.replication_factor_b)
                        y_f_correction = - 1.0 * parameters.replication_factor_b; 
                    if (guestbeads[first_beadid].z_f > 1.0 * parameters.replication_factor_c)
                        z_f_correction = - 1.0 * parameters.replication_factor_c; 
                    if (guestbeads[first_beadid].x_f < 0.0)
                        x_f_correction = 1.0 * parameters.replication_factor_a;
                    if (guestbeads[first_beadid].y_f < 0.0)
                        y_f_correction = 1.0 * parameters.replication_factor_b;
                    if (guestbeads[first_beadid].z_f < 0.0)
                        z_f_correction = 1.0 * parameters.replication_factor_c;

                    for (int b = 0; b < guestmoleculeinfo[which_type].nbeads; b++) {
                        int beadid = guestmolecules[idx_guestmolecule].beadID[b];
                        // adjust fractional coordinates, then later Cartesian coords
                        guestbeads[beadid].x_f += x_f_correction;
                        guestbeads[beadid].y_f += y_f_correction;
                        guestbeads[beadid].z_f += z_f_correction;
                        double x, y, z;
                        FractionalToCartesian(parameters.t_matrix, 
                                             guestbeads[beadid].x_f, guestbeads[beadid].y_f, guestbeads[beadid].z_f, 
                                             x, y, z); // TODO pass directly to here
                        guestbeads[beadid].x = x;
                        guestbeads[beadid].y = y;
                        guestbeads[beadid].z = z;
                    }
                    assert(OutsideUnitCell(guestbeads[first_beadid].x_f, 
                                           guestbeads[first_beadid].y_f, 
                                           guestbeads[first_beadid].z_f, 
                                           parameters) == false);

                    // get energy at new position
                    double E_gf_new = GuestFrameworkEnergy(idx_guestmolecule, 
                                guestmoleculeinfo,
                                guestmolecules,
                                guestbeads,
                                grid_info,
                                energy_grids);  // framework-guest
                    double E_gg_new = GuestGuestEnergy(N_g_total,
                            idx_guestmolecule, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            parameters);  // guest-guest
                    double E_new = E_gf_new + E_gg_new;
                    
                    // accept if, loosely, energeticall favorable
                    if (rand_for_acceptance < exp(-(E_new - E_old)/parameters.T)) {
                        stats.N_moves += 1; 
                        E_gg_this_cycle += E_gg_new - E_gg_old; E_gf_this_cycle += E_gf_new - E_gf_old;
                        if (E_new > 1e6)
                            std::cout << "Move accepted with huge energy" << std::endl;
                        // already overwrote coords with new coords, no need to update coords in this case
                    }
                    else {
                        // change new, moved cords back to old coords
                        for (int b = 0; b < guestmoleculeinfo[which_type].nbeads; b++) {
                            int beadid = guestmolecules[idx_guestmolecule].beadID[b];
                            guestbeads[beadid] = old_beads[b];
                        }
                    }
                } // end if N_g == 0
            } // end translation

            //
            //  PARTICLE IDENTITY CHANGE FOR DUAL COMPONENT
            //
            else if (whichMCmove < parameters.p_exchange + parameters.p_move + parameters.p_identity_change) {
                // if there are paricles of this type in the system
                if (N_g[which_type] > 0) { 
                    stats.N_ID_swap_trials += 1;

                    // pick which particles to change identity of
                    decltype(uniformint.param()) new_range (0, N_g[which_type] - 1); // set new range for rng
                    uniformint.param(new_range);
                    int idx_type = uniformint(generator); // which of this component?
                    int idx_guestmolecule = guestmolecule_index_list[which_type][idx_type]; // global index of this component
                    assert(guestmolecules[idx_guestmolecule].type == which_type);

                    // compute energy of guest with its current identity
                    double E_gf_old = GuestFrameworkEnergy(idx_guestmolecule, 
                                guestmoleculeinfo,
                                guestmolecules,
                                guestbeads,
                                grid_info,
                                energy_grids);  // framework-guest
                    double E_gg_old = GuestGuestEnergy(N_g_total,
                            idx_guestmolecule, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            parameters);  // guest-guest
                    double E_old = E_gf_old + E_gg_old;
                    
                    // store old guest and beads
                    GuestMolecule oldguestmolecule = guestmolecules[idx_guestmolecule];
                    GuestBead oldguestbeads[2];
                    for (int b = 0; b < guestmoleculeinfo[which_type].nbeads; b ++)
                        oldguestbeads[b] = guestbeads[oldguestmolecule.beadID[b]];

                    // in case there are two beads to start, declare which we start with.
                    int whichbead;
                   
                    // for ease, store number of beads in both and new_type 
                    int new_type = (which_type == 0) ? 1 : 0;
                    int nbeads_this = guestmoleculeinfo[which_type].nbeads;
                    int nbeads_new = guestmoleculeinfo[new_type].nbeads;
                    
                    // switch identity of guest
                    guestmolecules[idx_guestmolecule].type = new_type;

                    // single bead to single bead
                    if ((nbeads_this == 1) & (nbeads_new == 1)) {
                        // change type of bead
                        guestbeads[oldguestmolecule.beadID[0]].type = guestmoleculeinfo[new_type].beadtypes[0];
                    }

                    // single bead to double bead
                    if ((nbeads_this == 1) & (nbeads_new == 2)) {
                        // change identity of old (single) bead to the first bead of the new guest with two beads
                        guestbeads[oldguestmolecule.beadID[0]].type = guestmoleculeinfo[new_type].beadtypes[0];

                        // create and place the second bead
                        double theta = 2 * M_PI * uniform01(generator);
                        double phi = acos(2 * uniform01(generator) - 1);

                        // insert second bead at Cartesian coords
                        guestbeads[N_beads].x = guestbeads[guestmolecules[idx_guestmolecule].beadID[0]].x + 
                                                guestmoleculeinfo[new_type].bondlength * sin(phi) * cos(theta);
                        guestbeads[N_beads].y = guestbeads[guestmolecules[idx_guestmolecule].beadID[0]].y +
                                                guestmoleculeinfo[new_type].bondlength * sin(phi) * sin(theta);
                        guestbeads[N_beads].z = guestbeads[guestmolecules[idx_guestmolecule].beadID[0]].z + 
                                                guestmoleculeinfo[new_type].bondlength * cos(phi);

                        // convert to fractional 
                        double x_f, y_f, z_f;
                        CartesianToFractional(parameters.inv_t_matrix, 
                                              x_f, y_f, z_f, 
                                              guestbeads[N_beads].x, guestbeads[N_beads].y, guestbeads[N_beads].z);
                        // Remember, second bead allowed outside of bounding box...
                        // assign coords
                        guestbeads[N_beads].x_f = x_f;
                        guestbeads[N_beads].y_f = y_f;
                        guestbeads[N_beads].z_f = z_f;
                        // assign type and guestmolculeID
                        guestbeads[N_beads].type = guestmoleculeinfo[new_type].beadtypes[1];
                        guestbeads[N_beads].guestmoleculeID = idx_guestmolecule;

                        // make guestmolecule point to this new bead
                        guestmolecules[idx_guestmolecule].beadID[1] = N_beads;
                        // update bead count
                        N_beads += 1;
                    }

                    // double bead to single bead
                    if ((nbeads_this == 2) & (nbeads_new == 1)) {
                        // choose which bead to use for guest, randomly
                        whichbead = beadpicker(generator);
                        printf("which bead %d\n", whichbead);  // for testing

                        // change type of bead
                        guestbeads[oldguestmolecule.beadID[whichbead]].type = guestmoleculeinfo[new_type].beadtypes[0];
                        // make guest molecule point to this bead (energy computations will only see the first)
                        guestmolecules[idx_guestmolecule].beadID[0] = oldguestmolecule.beadID[whichbead];  // energy computations will ignore other bead...
                        // update bead count
                        N_beads -= 1;
                    }
                    
                    // double bead to double bead (deal with different lengths)
                    if ((nbeads_this == 2) & (nbeads_new == 2)) {
                        // choose which bead to hv in common
                        whichbead = beadpicker(generator);
                        int pivotbeadid = oldguestmolecule.beadID[whichbead];
                        // change type of bead, keep position
                        guestbeads[pivotbeadid].type = guestmoleculeinfo[new_type].beadtypes[0];
                        // insert other bead along same vector, just shrink length
                        //
                        int otherbead = (whichbead == 0) ? 1 : 0;
                        int otherbeadid = oldguestmolecule.beadID[otherbead];
                        // change type of other bead, but will move...
                        guestbeads[otherbeadid].type = guestmoleculeinfo[new_type].beadtypes[1];
                        assert(otherbead != whichbead);
                        
                        VectorR3 bondvector;  // get bond vector
                        bondvector.x = guestbeads[otherbeadid].x - guestbeads[pivotbeadid].x;
                        bondvector.y = guestbeads[otherbeadid].y - guestbeads[pivotbeadid].y;
                        bondvector.z = guestbeads[otherbeadid].z - guestbeads[pivotbeadid].z;

                        // move second bead along bond vector, closer/further to/from pivot bead
                        double bondlength_ratio = guestmoleculeinfo[new_type].bondlength / guestmoleculeinfo[which_type].bondlength;
                        guestbeads[otherbeadid].x = guestbeads[pivotbeadid].x + bondvector.x * bondlength_ratio;
                        guestbeads[otherbeadid].y = guestbeads[pivotbeadid].y + bondvector.y * bondlength_ratio;
                        guestbeads[otherbeadid].z = guestbeads[pivotbeadid].z + bondvector.z * bondlength_ratio;

                        double x_f, y_f, z_f;
                        CartesianToFractional(parameters.inv_t_matrix, 
                                             x_f, y_f, z_f, 
                                             guestbeads[otherbeadid].x, guestbeads[otherbeadid].y, guestbeads[otherbeadid].z);
                        guestbeads[otherbeadid].x_f = x_f;
                        guestbeads[otherbeadid].y_f = y_f;
                        guestbeads[otherbeadid].z_f = z_f;
                    }
                   
                    // OKAY, changed guestmolecules[idx_guestmolecule] type at this point... 
                    // if FIRST BEAD is moved outside of box, reflect ENTIRE guest to the other side
                    int first_beadid = guestmolecules[idx_guestmolecule].beadID[0];
                    double x_f_correction = 0.0; double y_f_correction = 0.0; double z_f_correction = 0.0;
                    if (guestbeads[first_beadid].x_f > 1.0 * parameters.replication_factor_a)
                        x_f_correction = - 1.0 * parameters.replication_factor_a;
                    if (guestbeads[first_beadid].y_f > 1.0 * parameters.replication_factor_b)
                        y_f_correction = - 1.0 * parameters.replication_factor_b;
                    if (guestbeads[first_beadid].z_f > 1.0 * parameters.replication_factor_c)
                        z_f_correction = - 1.0 * parameters.replication_factor_c;
                    if (guestbeads[first_beadid].x_f < 0.0)
                        x_f_correction = 1.0 * parameters.replication_factor_a;
                    if (guestbeads[first_beadid].y_f < 0.0)
                        y_f_correction = 1.0 * parameters.replication_factor_b;
                    if (guestbeads[first_beadid].z_f < 0.0)
                        z_f_correction = 1.0 * parameters.replication_factor_c;

                    for (int b = 0; b < guestmoleculeinfo[new_type].nbeads; b++) {
                        int beadid = guestmolecules[idx_guestmolecule].beadID[b];
                        // adjust fractional coordinates, then later Cartesian coords
                        guestbeads[beadid].x_f += x_f_correction;
                        guestbeads[beadid].y_f += y_f_correction;
                        guestbeads[beadid].z_f += z_f_correction;
                        double x, y, z;
                        FractionalToCartesian(parameters.t_matrix, 
                                             guestbeads[beadid].x_f, guestbeads[beadid].y_f, guestbeads[beadid].z_f, 
                                             x, y, z); // TODO pass directly to here
                        guestbeads[beadid].x = x;
                        guestbeads[beadid].y = y;
                        guestbeads[beadid].z = z;
                        assert(OutsideUnitCell(guestbeads[first_beadid].x_f, 
                                               guestbeads[first_beadid].y_f, 
                                               guestbeads[first_beadid].z_f, 
                                               parameters) == false);
                    }

                    // compute energy with the different identity
                    double E_gf_new = GuestFrameworkEnergy(idx_guestmolecule, 
                                guestmoleculeinfo,
                                guestmolecules,
                                guestbeads,
                                grid_info,
                                energy_grids);  // framework-guest
                    double E_gg_new = GuestGuestEnergy(N_g_total,
                            idx_guestmolecule, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            parameters);  // guest-guest
                    double E_new = E_gf_new + E_gg_new;
                    
                    // Accept move if, loosely, particle identity change was energetically favorable
                    double prob_acceptance_ID_swap = exp(-(E_new - E_old) / parameters.T) * parameters.fugacity[new_type] / parameters.fugacity[which_type] * static_cast<double>(N_g[which_type]) / (static_cast<double>( N_g[new_type]) + 1.0);
                    if (rand_for_acceptance < prob_acceptance_ID_swap) { 
                        stats.N_ID_swaps += 1;
                        
                        // keep track of energies
                        E_gg_this_cycle += E_gg_new - E_gg_old;
                        E_gf_this_cycle += E_gf_new - E_gf_old;
                        
                        // move last one of which_type to here, this one is deleted.
                        guestmolecule_index_list[which_type][idx_type] = guestmolecule_index_list[which_type][N_g[which_type] - 1]; 
                        // now this index is part of the other component ID
                        guestmolecule_index_list[new_type][N_g[new_type]] = idx_guestmolecule; 

                        // update particle numbers
                        N_g[which_type] -= 1; // one less of which type
                        N_g[new_type] += 1; // one more of new type
                        
                        // actually replace this deleted bead with last bead in array
                        if ((nbeads_this == 2) & (nbeads_new == 1)) { 
                            // replace this deleted guest's beads with the beads at the end of guestbeads
                            int otherbead = (whichbead == 0) ? 1 : 0;
                            printf("\tother bead%d\n", otherbead); // for test
                            int beadid_to_remove = oldguestmolecule.beadID[otherbead];
                            // move last bead to here. incremented N_beads =- 1 already
                            guestbeads[beadid_to_remove] = guestbeads[N_beads];
                            // change molecule ID of beads for this guest
                            int otherguestmoleculeid =  guestbeads[N_beads].guestmoleculeID;
                            // change beadID in this other guest molecule to beadID_to_remove
                            bool found = false;
                            for (int bi = 0; bi < guestmoleculeinfo[guestmolecules[otherguestmoleculeid].type].nbeads; bi++) {
                                if (guestmolecules[otherguestmoleculeid].beadID[bi] == N_beads) {
                                    guestmolecules[otherguestmoleculeid].beadID[bi] = beadid_to_remove;
                                    found = true;
                                }
                            }
                            assert(found == true); // TODO remove later, for checking
                        }
                    }  // end "if we accept this move"
                    else { 
                        // if didn't accept
                        // go back to old guest molecule position
                        guestmolecules[idx_guestmolecule] = oldguestmolecule;
                        // replace beads with old beads
                        for (int b = 0; b < guestmoleculeinfo[which_type].nbeads; b++)
                            guestbeads[oldguestmolecule.beadID[b]] = oldguestbeads[b];
                        
                        // ignore the last bead that we inserted by changing the number of beads
                        if ((nbeads_this == 1) & (nbeads_new == 2))
                            N_beads -= 1;
                        // change beads back to original value
                        if ((nbeads_this == 2) & (nbeads_new == 1))
                            N_beads += 1;
                    }  // end "if we didnt accept this move"
                } // end if N_g > 0
            } // end particle identity swap

            //
            //  MC trial: Regrow (move molecule to a different, random location, like a move mostly)
            //
            else {
                if (parameters.debugmode) 
                    printf("Regrow Trial.\n");
                stats.N_regrow_trials += 1;
                
                if (N_g[which_type] > 0) { 
                    // if there are paricles of this type in the system
                    // Randomly choose an adsorbate of which_type
                    decltype(uniformint.param()) new_range(0, N_g[which_type] - 1); // set new range for rng
                    uniformint.param(new_range);
                    int idx_regrow_type = uniformint(generator);  // corresponds to index of this type in guestmolecule_index
                    int idx_guestmolecule = guestmolecule_index_list[which_type][idx_regrow_type]; // global ID in guests
                    
                    // compute energy in current (old) position
                    double E_gf_old = GuestFrameworkEnergy(idx_guestmolecule, 
                                guestmoleculeinfo,
                                guestmolecules,
                                guestbeads,
                                grid_info,
                                energy_grids);  // framework-guest
                    double E_gg_old = GuestGuestEnergy(N_g_total,
                            idx_guestmolecule, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            parameters);  // guest-guest
                    double E_old = E_gf_old + E_gg_old;

                    GuestBead old_beads[2];  // for old coords of beads
                    // insert this guest at a new position
                    for (int b = 0; b < guestmoleculeinfo[which_type].nbeads; b++) {
                        int beadid = guestmolecules[idx_guestmolecule].beadID[b];
                        // store old coords of beads
                        old_beads[b] = guestbeads[beadid];
                        
                        // randomly insert first bead
                        if (b == 0) {
                            guestbeads[beadid].x_f = uniform01(generator) * parameters.replication_factor_a;
                            guestbeads[beadid].y_f = uniform01(generator) * parameters.replication_factor_b;
                            guestbeads[beadid].z_f = uniform01(generator) * parameters.replication_factor_c;
                            
                            double x, y, z; 
                            FractionalToCartesian(parameters.t_matrix, 
                                                  guestbeads[beadid].x_f, guestbeads[beadid].y_f, guestbeads[beadid].z_f, 
                                                  x, y, z);
                            guestbeads[beadid].x = x; 
                            guestbeads[beadid].y = y; 
                            guestbeads[beadid].z = z; 
                        }
                        // if 2nd bead, insert on sphere centered at first bead
                        if (b == 1) {
                            double theta = 2 * M_PI * uniform01(generator);
                            double phi = acos(2 * uniform01(generator) - 1);

                            // insert second bead at Cartesian coords
                            guestbeads[beadid].x = guestbeads[guestmolecules[idx_guestmolecule].beadID[0]].x + 
                                                    guestmoleculeinfo[which_type].bondlength * sin(phi) * cos(theta);
                            guestbeads[beadid].y = guestbeads[guestmolecules[idx_guestmolecule].beadID[0]].y +
                                                    guestmoleculeinfo[which_type].bondlength * sin(phi) * sin(theta);
                            guestbeads[beadid].z = guestbeads[guestmolecules[idx_guestmolecule].beadID[0]].z + 
                                                    guestmoleculeinfo[which_type].bondlength * cos(phi);

                            // convert to fractional to make sure we are inside bounding box
                            double x_f, y_f, z_f;
                            CartesianToFractional(parameters.inv_t_matrix, 
                                                  x_f, y_f, z_f, 
                                                  guestbeads[beadid].x, guestbeads[beadid].y, guestbeads[beadid].z);
                            
                            // SECOND bead is allowed outside simulation box. do not worry about PBC here!  assign fractional coords
                            guestbeads[beadid].x_f = x_f;
                            guestbeads[beadid].y_f = y_f;
                            guestbeads[beadid].z_f = z_f;
                        }  // end "if second bead..."
                    }  // end loop over beads

                    // get energy at new position
                    double E_gf_new = GuestFrameworkEnergy(idx_guestmolecule, 
                                guestmoleculeinfo,
                                guestmolecules,
                                guestbeads,
                                grid_info,
                                energy_grids);  // framework-guest
                    double E_gg_new = GuestGuestEnergy(N_g_total,
                            idx_guestmolecule, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            parameters);  // guest-guest
                    double E_new = E_gf_new + E_gg_new;
                    
                    // accept if, loosely, energeticall favorable
                    if (rand_for_acceptance < exp(-(E_new - E_old)/parameters.T)) {
                        stats.N_regrows += 1; 
                        E_gg_this_cycle += E_gg_new - E_gg_old; E_gf_this_cycle += E_gf_new - E_gf_old;
                        if (E_new > 1e6)
                            std::cout << "Move accepted with huge energy" << std::endl;
                        // already overwrote coords with new coords, no need to update coords in this case
                    }
                    else {
                        // change new, moved cords back to old coords
                        for (int b = 0; b < guestmoleculeinfo[which_type].nbeads; b++) {
                            int beadid = guestmolecules[idx_guestmolecule].beadID[b];
                            guestbeads[beadid] = old_beads[b];
                        }
                    }
                } // end if N_g == 0
            }

            // assert N_g < MAX_GUESTS
            if (N_g_total >= MAX_GUESTS - 2) {
                printf("N_g > MAX_GUESTS!\n");
                exit(EXIT_FAILURE);
            }

            //
            // Collect statistics
            //
            if ((cycle > parameters.numequilibriumtrials) & (cycle_counter % parameters.samplefrequency == 0)) {
                stats.N_samples += 1;
                stats.N_g_avg[0] += N_g[0]; 
                stats.N_g_avg[1] += N_g[1];
                stats.N_g2_avg[0] += N_g[0] * N_g[0]; 
                stats.N_g2_avg[1] += N_g[1] * N_g[1];
                stats.guest_guest_energy_avg += E_gg_this_cycle; // divide by N_samples later
                stats.framework_guest_energy_avg += E_gf_this_cycle; // divide by N_samples later
            }

            //
            // Print stuff if debug mode
            //
            if (parameters.debugmode) {
                // print guest molecule information
                printf("\n\n");
                printf("GUEST LIST OF %d guests:\n", N_g_total);
                for (int g = 0; g < N_g_total; g++) {
                    printf("\t%d. type %d. %d beads. Bead IDs:", g, guestmolecules[g].type, guestmoleculeinfo[guestmolecules[g].type].nbeads);
                    for (int b = 0; b < guestmoleculeinfo[guestmolecules[g].type].nbeads; b++)
                        printf(" %d", guestmolecules[g].beadID[b]);
                    printf("\n");
                }
                // print bead information
                printf("\nBEAD LIST OF %d beads.\n", N_beads);
                for (int b = 0; b < N_beads; b++) 
                    printf("\t%d. type %d, guestmolecule ID %d\n", b, guestbeads[b].type, guestbeads[b].guestmoleculeID);
                // print adsorbate list
                printf("GUEST LIST BY TYPE.\n");
                printf("\tType 0 guest IDs: ");
                for (int a = 0 ; a < N_g[0] ; a++)
                     printf(" %d", guestmolecule_index_list[0][a]);
                if (parameters.numadsorbates > 1) {
                    printf("\n\tType 1 guest IDs: ");
                    for (int a = 0 ; a < N_g[1] ; a++)
                         printf(" %d", guestmolecule_index_list[1][a]);
                }
                printf("\n\n");
            }
            
            //
            // Make assertions to check code
            //
            if (parameters.makeassertions) {
                double bondlength_tol = 0.000001;  // tolerance for bond length
                double coords_tol = 0.000001; // tolerance for coordinates differences
                int typecounts[2] = {0};
                int beadcount = 0;
                for (int g = 0; g < N_g_total; g++) {
                    // count molecules of each type
                    typecounts[guestmolecules[g].type] ++;
                    // count beads
                    beadcount += guestmoleculeinfo[guestmolecules[g].type].nbeads;
                    // assert first bead is inside unit cell
                    int first_beadid = guestmolecules[g].beadID[0];
                    assert(OutsideUnitCell(guestbeads[first_beadid].x_f, guestbeads[first_beadid].y_f, guestbeads[first_beadid].z_f, parameters) == false);
                    // assert bond lengths if applicable
                    if (guestmoleculeinfo[guestmolecules[g].type].nbeads == 2) {
                        int second_beadid = guestmolecules[g].beadID[1];
                        double dx = guestbeads[first_beadid].x - guestbeads[second_beadid].x;
                        double dy = guestbeads[first_beadid].y - guestbeads[second_beadid].y;
                        double dz = guestbeads[first_beadid].z - guestbeads[second_beadid].z;
                        double l = sqrt(dx*dx + dy*dy + dz*dz);
                        assert(l < guestmoleculeinfo[guestmolecules[g].type].bondlength + bondlength_tol);
                        assert(l > guestmoleculeinfo[guestmolecules[g].type].bondlength - bondlength_tol);
                    }
                    // assert beads in this guest pointing correctly to this guest molecule
                    for (int b = 0; b < guestmoleculeinfo[guestmolecules[g].type].nbeads; b++) {
                        int beadid = guestmolecules[g].beadID[b];
                        if (g != guestbeads[beadid].guestmoleculeID) {
                            printf("FAIL, guest bead %d is assigned to guest molecule %d\n", b, guestbeads[b].guestmoleculeID);
                            printf("\tbut bead IDs of guest molecule %d: %d and %d\n", g, guestmolecules[g].beadID[0], guestmolecules[g].beadID[1]);
                            exit(EXIT_FAILURE);
                        }
                    }
                }  // end loop over guests
                assert(beadcount == N_beads);
                assert(typecounts[0] == N_g[0]);
                assert(typecounts[1] == N_g[1]);
                assert(N_g[0] + N_g[1] == N_g_total);
                
                // assert that our index list is correct
                for (int i0 = 0; i0 < N_g[0]; i0 ++) {
                    int idx_guestmolecule = guestmolecule_index_list[0][i0]; 
                    assert(guestmolecules[idx_guestmolecule].type == 0);
                }
                for (int i1 = 0; i1 < N_g[1]; i1 ++) {
                    int idx_guestmolecule = guestmolecule_index_list[1][i1]; 
                    assert(guestmolecules[idx_guestmolecule].type == 1);
                }

                // assert guest beads fractional and cartesian coords match
                for (int b = 0; b < N_beads; b ++) {
                    double x, y, z;
                    FractionalToCartesian(parameters.t_matrix,
                                         guestbeads[b].x_f, guestbeads[b].y_f, guestbeads[b].z_f,
                                         x, y, z);

                    if (
                    (guestbeads[b].x > x + coords_tol) |
                    (guestbeads[b].y > y + coords_tol) |
                    (guestbeads[b].z > z + coords_tol) |
                    (guestbeads[b].x < x - coords_tol) |
                    (guestbeads[b].y < y - coords_tol) |
                    (guestbeads[b].z < z - coords_tol)
                    ) {
                        printf("FAIL on bead id %d, belongs to guest %d\n", b, guestbeads[b].guestmoleculeID);
                        printf("\tGuest molecule beadIDs: %d and %d", guestmolecules[guestbeads[b].guestmoleculeID].beadID[0], guestmolecules[guestbeads[b].guestmoleculeID].beadID[1]);
                        printf("\nGuest type %d, index %d, total guests %d\n", guestmolecules[guestbeads[b].guestmoleculeID].type, guestbeads[b].guestmoleculeID, N_g_total);
                        printf(" (x_f, y_f, z_f) = (%f, %f, %f)\n" , guestbeads[b].x_f, guestbeads[b].y_f, guestbeads[b].z_f);
                        printf(" (x, y, z) = (%f, %f, %f)\n" , guestbeads[b].x, guestbeads[b].y, guestbeads[b].z);
                        printf(" (x, y, z) = (%f, %f, %f)\n" , x, y, z);
                        exit(EXIT_FAILURE);
                    }
               }
            }  // end make assertions
            
            //
            // Write adsorbate positions to file (optional)
            //
            if (parameters.writeadsorbatepositions) {
                if ((cycle > parameters.numequilibriumtrials) & (cycle_counter % parameters.writepositionfrequency == 0)) {
                    N_snapshots ++;
                    WriteGuestPostionsToFile(adsorbatepositionfile, 
                              N_g_total,
                              guestmoleculeinfo,
                              guestmolecules,
                              guestbeads,
                              parameters); 
                    if (N_snapshots > parameters.num_snapshots) {
                        printf("Reached %d snapshots, exiting.\n", N_snapshots);
                        fprintf(outputfile, "\nWrote %d adsorbate snapshot positions in xyz file every %d MC moves.\n", 
                                        N_snapshots, parameters.writepositionfrequency);
                        std::exit(EXIT_SUCCESS);
                    }
                }
            }

        }  // end inner cycle loop
//
    }  // end outer cycle loop
    
    // take avg
    stats.guest_guest_energy_avg = 1.0 * stats.guest_guest_energy_avg / stats.N_samples;
    stats.framework_guest_energy_avg = 1.0 * stats.framework_guest_energy_avg / stats.N_samples;

    double sim_time = ReadTimer() - start_of_sim_time;
    fprintf(outputfile, "\nEnergy checks\n");
    // check energy calcs
    double E_gg_system = TotalGuestGuest(N_g_total,
                        guestmoleculeinfo,
                        guestmolecules,
                        guestbeads,
                        parameters);
    double E_gf_system = TotalGuestFrameworkEnergy(N_g_total, 
                            guestmoleculeinfo,
                            guestmolecules,
                            guestbeads,
                            grid_info,
                            energy_grids);
    fprintf(outputfile, "    E_gg total calc'ed at end: %f K\n", E_gg_system);
    fprintf(outputfile, "    E_gg from adding dE's throughout simulation: %f K\n", E_gg_this_cycle);
    fprintf(outputfile, "    E_gf total calc'ed at end: %f K\n", E_gf_system);
    fprintf(outputfile, "    E_gf from adding dE's throughout simulation: %f K\n", E_gf_this_cycle);
    
    // write stats
    fprintf(outputfile, "\nGCMC Statistics\n");
    fprintf(outputfile, "    Simulation time: %f min\n", sim_time / 60.0);
    fprintf(outputfile, "    Insertions: %d / %d (%f %% accepted)\n", stats.N_insertions, stats.N_insertion_trials, 100.0 * stats.N_insertions / stats.N_insertion_trials);
    fprintf(outputfile, "    Deletions: %d / %d (%f %% accepted)\n", stats.N_deletions, stats.N_deletion_trials, 100.0 * stats.N_deletions / stats.N_deletion_trials);
    if (parameters.p_move > 0.0)
        fprintf(outputfile, "    Moves: %d / %d (%f %% accepted)\n", stats.N_moves, stats.N_move_trials, 100.0 * stats.N_moves / stats.N_move_trials);
    if (parameters.p_identity_change > 0.0)
        fprintf(outputfile, "    Identity changes: %d / %d (%f %% accepted)\n", stats.N_ID_swaps, stats.N_ID_swap_trials, 100.0 * stats.N_ID_swaps / stats.N_ID_swap_trials);
    if (parameters.p_regrow > 0.0)
        fprintf(outputfile, "    Regrows: %d / %d (%f %% accepted)\n", stats.N_regrows, stats.N_regrow_trials, 100.0 * stats.N_regrows / stats.N_regrow_trials);
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


    fclose(outputfile); 
}
