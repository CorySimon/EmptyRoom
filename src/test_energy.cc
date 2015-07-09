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
        printf("Bead %d: energy = %f\n", b, BeadFrameworkEnergy(WrapToInterval(guestbeads[beadid].x_f, 1.0),
                                            WrapToInterval(guestbeads[beadid].y_f, 1.0),
                                            WrapToInterval(guestbeads[beadid].z_f, 1.0),
                                            grid_info, energy_grids[beadtype]));
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
    // Write settings to outputfile
    //
    FILE * outputfile;
    outputfile = fopen("test_energy.txt", "w");
    WriteSettingsToOutputfile(outputfile, parameters, framework, forcefield, grid_info, guestmoleculeinfo);
    if (parameters.verbose) printf("Wrote info to outputfile\n");

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
    
    // to test energy.
    guestmolecules[0].type = 0;
    guestmolecules[0].beadID[0] = 0;
    guestbeads[0].x_f = 0.0000001;
    guestbeads[0].y_f = 0.0000001;
    guestbeads[0].z_f = 0.0000001;
    double x, y, z; 
    FractionalToCartesian(parameters.t_matrix, 
                          guestbeads[0].x_f, guestbeads[0].y_f, guestbeads[0].z_f, 
                          x, y, z);
    guestbeads[0].x = x; 
    guestbeads[0].y = y; 
    guestbeads[0].z = z; 
    // define which guest this bead belongs to
    guestbeads[0].guestmoleculeID = 0;
    // define bead type for energy computations
    guestbeads[0].type = guestmoleculeinfo[0].beadtypes[0];

    // add second bead if not LJ sphere
    // generate Cartesian coords of second bead
    guestbeads[1].x = guestbeads[0].x + guestmoleculeinfo[0].bondlength; 
    guestbeads[1].y = guestbeads[0].y;
    guestbeads[1].z = guestbeads[0].z;

    // convert to fractional 
    double x_f, y_f, z_f;
    CartesianToFractional(parameters.inv_t_matrix, 
                          x_f, y_f, z_f, 
                          guestbeads[1].x, guestbeads[1].y, guestbeads[1].z);
    guestbeads[1].x_f = x_f;
    guestbeads[1].y_f = y_f;
    guestbeads[1].z_f = z_f;  // TODO pass this directly to cartesiantofrac
    
    // SECOND BEAD IS ALLOWED OUTSIDE OF THE UNIT CELL! (not the first though) 
     
    // define which guest this bead belongs to
    guestbeads[1].guestmoleculeID = 0;
    // define bead type for energy computations
    guestbeads[1].type = guestmoleculeinfo[0].beadtypes[1];
    
    // add this bead to guests
    guestmolecules[0].beadID[1] = 1;
    double E_gf = GuestFrameworkEnergy(0, 
                guestmoleculeinfo,
                guestmolecules,
                guestbeads,
                grid_info,
                energy_grids);  // framework-guest
    
    for (int i=0; i<2; i++)
        printf("Bead %d: [%f, %f, %f] fractional.\n", i, guestbeads[i].x_f, guestbeads[i].y_f, guestbeads[i].z_f);
    printf("E_gf = %f\n", E_gf);
    
}
