//
//  GCMC code (mu, V, T) simulations.
//
#include <chrono>
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include <limits>
#include<set>
#include<map>
#include <complex>
#include <cstring>
#include<math.h>
#include<random>
#include<assert.h>
#include<cstdlib> // for "exit"
#include<sstream> // string stream
#include<vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
namespace boo = boost::numeric::ublas;
#include "adsorbate.h"
#include "datatypes.h"
#include "Framework.h"
#include "Forcefield.h"
#include <sys/time.h>
#include "readsettings.h"
#include "Forcefield.h"
#include "write_settings_to_outputfile_gcmc.h"
#include <sys/time.h>
//#include "pocketblocking.h"
#define min_r .0000000000000001  // don't want to divide by zero...

/* 
 * Global Variables for ease. BE CAREFUL!
*/
boo::matrix<double> t_matrix(3, 3);  // fractional to Cartesian coords
boo::matrix<double> inv_t_matrix(3, 3);  // Cartesian to fractional

// storing force field parameters
boo::matrix<double> epsilon_matrix;
boo::matrix<double> sigma_squared_matrix;

// adsorbate names
std::vector<std::string> adsorbate;
// fugacity of adsorbate species
std::vector<double> fugacity;  // fugacity of adsorbates
// map of unique beads in adsorbates and bead type integers
std::map<std::string, int> beadlabel_to_int;
std::map<int, std::string> int_to_beadlabel;
std::vector<double> adsorbateMW;

// templates of adsorbate molecules
std::vector<Adsorbate> adsorbatetemplates;
// vector of adsorbate molecules for simulation
std::vector<Adsorbate> adsorbates;

// replication factor for unit cell
std::vector<int> uc_reps(3);

double r_cutoff_squared;  // (A), Lennard Jones cutoff radius, squared

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

double WrapToInterval(double x, double z) {
    // for applying periodic bc's
    return x - z * floor(x / z);
}

double GuestGuestEnergy(std::vector<Adsorbate> & adsorbates,
                        int adsorbateid) { 
    // Compute potential energy of adsorbate molecule adsorbateid in the vector adsorbates
    // Energy contribution from other guests.
    double E_gg = 0.0; // initiate
    
    // for each bead in this adsorbate molecule
    for (int b_this = 0; b_this < adsorbates[adsorbateid].nbeads; b_this++) {
        // get fractional coord of this bead in this adsorbate
        boo::matrix_column<boo::matrix<double> > xf_this_bead(adsorbates[adsorbateid].bead_xyz_f, b_this);
        int this_bead_type = adsorbates[adsorbateid].beadtypes[b_this];
        // for each other guest molecule 
        for (int k = 0 ; k < adsorbates.size(); k++) {
            // do not include self interation!
            if (k == adsorbateid) 
                continue; 
            // get energy contribution from each bead in other adsorbate molecule
            for (int b_other = 0; b_other < adsorbates[k].nbeads; b_other++) {
                int other_bead_type = adsorbates[k].beadtypes[b_other];
                boo::matrix_column<boo::matrix<double> > xf_other_bead(adsorbates[k].bead_xyz_f, b_other);

                // distance in fractional coords
                boo::vector<double> dx_f = xf_this_bead - xf_other_bead;

                // take nearest image
                for (int i_ = 0; i_ < 3; i_++) {
                    dx_f[i_] = dx_f[i_] - uc_reps[i_] * round(dx_f[i_] / uc_reps[i_]);
                    // assert within bounds #todo remove later
                    assert(dx_f[i_] < 0.5 * uc_reps[i_]);
                    assert(dx_f[i_] > -0.5 * uc_reps[i_]);
                }

                // distance in Cartesian
                boo::vector<double> dx = prod(t_matrix, dx_f);
                
                // Compute LJ potential
                double r2 = inner_prod(dx, dx);
                if (r2 < min_r) 
                    return 100000000000000.0; // overwrite if too small
                if (r2 < r_cutoff_squared) {
                    double sigma_over_r_sixth = pow(sigma_squared_matrix(this_bead_type, other_bead_type) / r2, 3.0); 
                    E_gg += 4.0 * epsilon_matrix(this_bead_type, other_bead_type) * sigma_over_r_sixth * (sigma_over_r_sixth - 1.0);
                }
            }  // end loop over BEADS of other guest molecule
        }  // end loop over other guest molecules
    }  // end loop over beads of THIS guest molecule
    return E_gg;
}

double TotalGuestGuestEnergy(std::vector<Adsorbate> & adsorbates) {
    // Calculate total system guest-guest energy
    double E_gg = 0.0;
    for (int i = 0; i < adsorbates.size(); i ++) 
        E_gg += GuestGuestEnergy(adsorbates, i);
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

double GuestFrameworkEnergy(Adsorbate adsorbate,
                            GridInfo grid_info,
                            double ** energy_grids) {
    // Compute energy of adsorbate with framework
    double E_gf = 0.0;
    for (int b = 0; b < adsorbate.nbeads; b++) {  
        // for each bead that makes up this guest molecule
        E_gf += BeadFrameworkEnergy(WrapToInterval(adsorbate.bead_xyz_f(0, b), 1.0),
                                    WrapToInterval(adsorbate.bead_xyz_f(1, b), 1.0),
                                    WrapToInterval(adsorbate.bead_xyz_f(2, b), 1.0),
                                    grid_info, energy_grids[adsorbate.beadtypes[b]]);
    }
    return E_gf;
}

double TotalGuestFrameworkEnergy(std::vector<Adsorbate> & adsorbates,
                                 GridInfo grid_info,
                                 double ** energy_grids) {
    // Calculate total system framework-guest energy
    double E_gf = 0.0;
    for (int i = 0; i < adsorbates.size(); i ++)
        E_gf += GuestFrameworkEnergy(adsorbates[i], grid_info, energy_grids);
    return E_gf;
}

////
//// Write guest positions to file
////
//void WriteGuestPostionsToFile(FILE * positionfile, 
//                              int N_g_total,
//                              GuestMoleculeInfo * guestmoleculeinfo,
//                              GuestMolecule * guestmolecules,
//                              GuestBead * guestbeads,
//                              GCMCParameters parameters) {
//    // for each guest molecule...
//    for (int i = 0; i < N_g_total ; i++) {
//        // for each bead in this guest molecule...
//        for (int b = 0; b < guestmoleculeinfo[guestmolecules[i].type].nbeads; b++) {
//            int beadid = guestmolecules[i].beadID[b];  // ID of this bead
//            int beadtype = guestbeads[beadid].type;
//            std::string beadlabel = (guestmoleculeinfo[guestmolecules[i].type].beadtypes[0] == beadtype) ? 
//                                    guestmoleculeinfo[guestmolecules[i].type].beadlabels[0] : 
//                                    guestmoleculeinfo[guestmolecules[i].type].beadlabels[1]; 
//            fprintf(positionfile, "%s %f %f %f\n", 
//                    beadlabel.c_str(),
//                    guestbeads[beadid].x, guestbeads[beadid].y, guestbeads[beadid].z);
//        }
//    }
//}
//

int main(int argc, char *argv[])
{
    if (! ((argc == 4) | (argc == 6) | (argc == 8))) {
        printf("Run as ./gcmc $structure $adsorbate0 $fugacity0(Pa) $adsorbate1 $fugactiy1(Pa) $adsorbate2 $fugactiy2(Pa)\nAdsorbate1,2 stuff is optional\n");
        exit(EXIT_FAILURE);
    }
    
    GCMCParameters parameters;
    // read arguments to get adsorbate and fugacity
    parameters.frameworkname = argv[1];

    if (argc == 4)
        parameters.numadsorbates = 1;
    if (argc == 6)
        parameters.numadsorbates = 2;
    if (argc == 8)
        parameters.numadsorbates = 3;
    for (int a = 0; a < parameters.numadsorbates; a++) {
        adsorbate.push_back(argv[2 + 2 * a]);
        fugacity.push_back(atof(argv[3 + 2 * a]));
    }

    ReadSimulationInputFile(parameters);
    if (parameters.verbose) 
        printf("Read simulation.input\n");
    if ((parameters.numadsorbates == 1) & (parameters.p_identity_change > 0.0)) {
        printf("Only 1 adsorbate and ID swap probability > 1... Make it zero.\n");
        exit(EXIT_FAILURE);
    }
    r_cutoff_squared = parameters.r_cutoff_squared;  // make global var
    
    // uc needs to be at least twice the cutoff radius for only methane within r_c to be within cutoff
    uc_reps = ReadUnitCellReplicationFile(parameters.frameworkname, "twice");
    if (parameters.verbose) 
        printf("Read .uc replication file\n");
    
    for (int i = 0; i < parameters.numadsorbates; i++) 
        adsorbateMW.push_back(GetAdsorbateMW(adsorbate[i]));
    
    //
    // Now that we know number of adsorbates, initialize guest count of each type
    //
    std::vector<int> N_g(parameters.numadsorbates, 0);

    //
    // Construct forcefield and framework objects
    //
    Forcefield forcefield(parameters.forcefieldname);
    if (parameters.verbose) 
        printf("Constructed Forcefield object\n");

    Framework framework(parameters.frameworkname);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            t_matrix(i, j) = framework.t_matrix[i][j];
            inv_t_matrix(i, j) = framework.inv_t_matrix[i][j];
        }
    }
    parameters.N_framework_atoms = framework.noatoms;
    if (parameters.verbose) 
        printf("Constructed Framework object\n");

    //
    // Read adsorbate info
    // Get list of unique beads and mapping from label to int
    //
    beadlabel_to_int = GetBeadMap(adsorbate, false);
    // create reverse map
    for (std::map<std::string, int>::iterator it=beadlabel_to_int.begin(); it!=beadlabel_to_int.end(); ++it)
        int_to_beadlabel[it->second] = it->first;
    adsorbatetemplates = GetAdsorbateTemplates(adsorbate, beadlabel_to_int, false);
    int nuniquebeads = beadlabel_to_int.size();
    if (parameters.verbose) {
        printf("Adsorbate templates read. %d unique beads.\n", nuniquebeads);
        for (int i = 0; i < adsorbate.size(); i++) {
            adsorbatetemplates[i].print_info();
        }
    }
    
    //
    // Get adsorbate LJ params for each bead 
    // Store in epsilon and sigma_squared matrices
    //
    epsilon_matrix.resize(nuniquebeads, nuniquebeads);  // pre-allocate size
    sigma_squared_matrix.resize(nuniquebeads, nuniquebeads);
    for (int bx = 0; bx < nuniquebeads; bx++) {
        std::vector<double> eps_sig_x = GrabGuestForceFieldParams(forcefield, int_to_beadlabel[bx]); 
        for (int by = 0; by < nuniquebeads; by++) {
            std::vector<double> eps_sig_y = GrabGuestForceFieldParams(forcefield, int_to_beadlabel[by]);
            epsilon_matrix(bx, by) = sqrt(eps_sig_x[0] * eps_sig_y[0]); 
            double sigma_mixed = (eps_sig_x[1] + eps_sig_y[1]) / 2.0;
            sigma_squared_matrix(bx, by) = sigma_mixed * sigma_mixed;
        }
    }
    if (parameters.verbose) 
        printf("Fetched adsorbate FF parameters\n");

    //
    // Import energy grid for each bead
    // + pocket blocking if enabled.
    //
    GridInfo grid_info; // for storing grid info
    double ** energy_grids = (double **) malloc(nuniquebeads * sizeof(double *));
    bool pocket_block_verbose_mode = false; // writes grid before and after
    for (int n_c = 0; n_c < nuniquebeads; n_c ++) {
        int N_x_temp, N_y_temp, N_z_temp;
        if (parameters.verbose) 
            printf("Importing energy grid %d\n", n_c);
        char gridfilename[512];
        sprintf(gridfilename, "data/grids/%s_%s_%s.txt", framework.name.c_str(), int_to_beadlabel[n_c].c_str(), forcefield.name.c_str());

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
            printf("energy grid %d for bead %s imported successfully.\n", n_c, int_to_beadlabel[n_c].c_str());
        
//        //
//        // Flood fill/ pocket blocking, if enabled
//        //
        if (parameters.pocketblocking) {
            printf("Pocket blocking not avail yet");
            exit(EXIT_FAILURE);
        }
//            double time_before = ReadTimer();
//            if (pocket_block_verbose_mode) {
//                printf("Pocket blocking beginning. Write a cube for before\n");
//                if (n_c == 0)
//                    WriteCube("before_blocking_0", framework, parameters, energy_grids[0], grid_info); //just so you can see the grid before ...
//            }
//
//            // pass energy grid to FindAndBlockPockets, which will do the job.
//            grid_info.numpockets[n_c] = FindAndBlockPockets(energy_grids[n_c], grid_info, parameters.T, parameters);
//            
//            if (pocket_block_verbose_mode) {
//                printf("Pocket blocking finished. Write a cube for after\n");
//                double time_after = ReadTimer();
//                printf("Time spent to find and block pockets: %f s\n", time_after - time_before);
//                if (n_c == 0)
//                    WriteCube("after_blocking_0", framework, parameters, energy_grids[0], grid_info); // ... and after pocket blocking
//            }
//        }
    }
    
    //
    // Initialize stats and guests array (includes both components via "type" attribute)
    //
    GCMCStats stats;
    InitializeGCMCStats(stats);
    double volume = framework.volume_unitcell * uc_reps[0] * uc_reps[1] * uc_reps[2]; // A ^ 3
//    outputfile << "\tSize of guests = " << MAX_GUESTS * sizeof(particle_g) / (1024.0 * 1024.0) << " MB\n";
    
    //
    // Set up random number generators
    //
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);  // Mersenne Twister algo
    std::uniform_real_distribution<double> uniform01(0.0, 1.0); // uniformly distributed real no in [0,1]
    std::uniform_int_distribution<int> uniformint(0, 10); // initialize with bogus range, will change later.
    std::uniform_int_distribution<int> type_generator(0, parameters.numadsorbates - 1); // for picking a component
    std::normal_distribution<double> std_normal_distn(0.0, 1.0);  // std normal for rotations: mean 0, std 1.0
    if (parameters.verbose) 
        printf("Random number generators initialized...\n");

    //
    // Write settings to outputfile
    //
    FILE * outputfile;
    char outputfilename[1024];
    
    // use argv instead of fugacity to make querying files easy
    if (parameters.numadsorbates == 1)
        sprintf(outputfilename, "output_files/%s_%s_%sPa_%gK_gcmc.out", 
            parameters.frameworkname.c_str(), adsorbate[0].c_str(), argv[3], parameters.T);
    if (parameters.numadsorbates == 2)
        sprintf(outputfilename, "output_files/%s_%s_%sPa_%s_%sPa_%gK_gcmc.out", 
            parameters.frameworkname.c_str(), 
            adsorbate[0].c_str(), argv[3],
            adsorbate[1].c_str(), argv[5],
            parameters.T);
    if (parameters.numadsorbates == 3)
        sprintf(outputfilename, "output_files/%s_%s_%sPa_%s_%sPa_%s_%sPa_%gK_gcmc.out", 
            parameters.frameworkname.c_str(), 
            adsorbate[0].c_str(), argv[3],
            adsorbate[1].c_str(), argv[5],
            adsorbate[2].c_str(), argv[7],
            parameters.T);
    outputfile = fopen(outputfilename, "w");
    WriteSettingsToOutputfile(outputfile,
        parameters,
        framework,
        forcefield,
        grid_info,
        uc_reps,
        t_matrix,
        inv_t_matrix,
        adsorbate,
        adsorbateMW,
        int_to_beadlabel,
        adsorbatetemplates,
        epsilon_matrix,
        sigma_squared_matrix);
    if (parameters.verbose) 
        printf("Wrote info to outputfile\n");

    //
    // Initialize adsorbate positions file to dump adsorbate positions
    // Only if option is enabled via parameters.writeadsorbatepositions bool
    //
    FILE * adsorbatepositionfile;
    int N_snapshots = 0;
    if (parameters.writeadsorbatepositions) {
        char positionfilename[1024];
        sprintf(positionfilename, "outputfiles/adsorbate_positions_%s.xyz", framework.name.c_str());
        adsorbatepositionfile = fopen(positionfilename, "w");
    }

    //
    //  Run GCMC
    //
    double start_of_sim_time = ReadTimer();
    if (parameters.verbose) 
        printf("Starting simulation...\n");
    double E_gg_this_cycle = 0.0; 
    double E_gf_this_cycle = 0.0; // assumes that we start with zero particles! (way it is set now...)
    int cycle_counter = 0;

    for (int cycle = 0; cycle < parameters.numtrials; cycle++) {
        // each cycle corresponds to a number of Markov chain moves defined by ninnercycles
        int ninnercycles = adsorbates.size() > 20 ? adsorbates.size() : 20; // proportional to number of guests
        for (int inner_cycle = 0; inner_cycle < ninnercycles; inner_cycle ++) {
            cycle_counter += 1;
            
            // generate random numbers
            double whichMCmove = uniform01(generator); // Insertion, deletion, translation, ID swap, or regrow?
            double rand_for_acceptance = uniform01(generator); // for testing acceptance
            int which_type = type_generator(generator); // select guest type 
            
            //
            //  MC trial: Insertion
            //
            if (whichMCmove < parameters.p_exchange / 2.0) {
                stats.N_insertion_trials += 1;

                // add new guest of this type from the template
                Adsorbate inserted_adsorbate = adsorbatetemplates[which_type];
                // if >1 beads, perform a random rotation
                if (inserted_adsorbate.nbeads > 1)
                    PerformUniformRandomRotation(inserted_adsorbate, generator, std_normal_distn);

                // translate to these fractional coords
                boo::vector<double> xf(3);
                xf[0] =  uniform01(generator) * uc_reps[0];
                xf[1] =  uniform01(generator) * uc_reps[1];
                xf[2] =  uniform01(generator) * uc_reps[2];
                // =  these Cartesian
                boo::vector<double> x = boo::prod(t_matrix, xf);
                
                // translate adsorbate to these coords (also updates fractional coords)
                inserted_adsorbate.translate_by_Cartesian_vector(x, t_matrix, inv_t_matrix, uc_reps);

                // append this new adsorbate to end of adsorabtes vector
                adsorbates.push_back(inserted_adsorbate);

                // compute energy of this inserted particle
                double E_gf = GuestFrameworkEnergy(inserted_adsorbate, grid_info, energy_grids);
                double E_gg = GuestGuestEnergy(adsorbates, adsorbates.size() - 1);
                double E_insertion = E_gf + E_gg;
                
                // for debug mode, print stuff off
                if (parameters.debugmode) {
                    printf("INSERTION PROPOSAL of adsorbate type %d\n", which_type);
                    printf("\n\tTranslate proposal at xf=%f,yf=%f,zf=%f\n", 
                                xf[0], xf[1], xf[2]);
                    if (inserted_adsorbate.nbeads > 1) {
                        boo::matrix_column<boo::matrix<double> > x_bead1(inserted_adsorbate.bead_xyz, 0);
                        boo::matrix_column<boo::matrix<double> > x_bead2(inserted_adsorbate.bead_xyz, 1);
                        printf("\tBond length between bead 0 and 1 = %f\n", boo::norm_2(x_bead1 - x_bead2));
                    }
                    std::cout << "\tE_gg = " << E_gg << std::endl;
                    std::cout << "\tE_gf = " << E_gf << std::endl;
                }
                
                // accept, loosely, if energetically favorable
                double acceptance_insertion = fugacity[which_type] * volume / ((N_g[which_type] + 1) * 1.3806488e7 * parameters.T) * exp(-E_insertion / parameters.T);
                if (parameters.debugmode) 
                    std::cout << "\tAcceptance prob = " << acceptance_insertion << std::endl;
                if (rand_for_acceptance < acceptance_insertion) {  
                    if (parameters.debugmode) 
                        printf("\tInsertion ACCEPTED.\n");

                    stats.N_insertions += 1;
                    
                    // update molecule and bead count
                    N_g[which_type] += 1;
                    
                    // update system energies
                    E_gg_this_cycle += E_gg; 
                    E_gf_this_cycle += E_gf;
                }
                else {
                    // remove adsorbate from vector
                    adsorbates.pop_back();
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
                    int which_of_this_adsorbate_type = uniformint(generator);  // which guest to propose to delete of this type

                    // loop through adsorbates vector until we find the which_of_this_adsorbate_type'th guest of this type
                    int count_this_type = 0;
                    int idx_delete;  // index of adsorbate in adsorbates we propose to delete
                    for (idx_delete = 0; idx_delete < adsorbates.size(); idx_delete++) {
                        if (adsorbates[idx_delete].type == which_type) {
                            if (count_this_type == which_of_this_adsorbate_type)
                                break;
                            count_this_type += 1;
                        }
                    }
                      
                    if (parameters.debugmode)
                        printf("\tProposal to delete guest molecule %d\n", idx_delete);
                    
                    // compute energy of this guest that we propose to delete
                    double E_gf = GuestFrameworkEnergy(adsorbates[idx_delete], grid_info, energy_grids);
                    double E_gg = GuestGuestEnergy(adsorbates, idx_delete);
                    double E_deletion = E_gf + E_gg;
                    
                    // accept deletion if, loosely, energetically favorable
                    double acceptance_del = (N_g[which_type] * 1.3806488e7 * parameters.T) / (fugacity[which_type] * volume) * exp(E_deletion / parameters.T);
                    if (rand_for_acceptance < acceptance_del) {
                        stats.N_deletions += 1;
                        
                        // deincrement count of guests of this species
                        N_g[which_type] -= 1;
                       
                        // erase this guest from the adsorbates vector
                        adsorbates.erase(adsorbates.begin() + idx_delete);

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
                    
                    // randomly select guest of this type to propose to delete
                    int which_of_this_adsorbate_type = uniformint(generator);  // which guest to propose to delete of this type

                    // loop through adsorbates vector until we find the which_of_this_adsorbate_type'th guest of this type
                    int count_this_type = 0;
                    int idx_move;  // index of adsorbate in adsorbates we propose to move
                    for (idx_move = 0; idx_move < adsorbates.size(); idx_move++) {
                        if (adsorbates[idx_move].type == which_type) {
                            if (count_this_type == which_of_this_adsorbate_type)
                                break;
                            count_this_type += 1;
                        }
                    }
                      
                    if (parameters.debugmode)
                        printf("\tProposal to move guest molecule %d\n", idx_move);
                    
                    // compute energy of this guest that we propose to move
                    double E_gf_old = GuestFrameworkEnergy(adsorbates[idx_move], grid_info, energy_grids);
                    double E_gg_old = GuestGuestEnergy(adsorbates, idx_move);
                    double E_old = E_gf_old + E_gg_old;

                    // store this adsorbate's attributes before it was moved
                    // This way, if move is rejected, we can just replace it
                    Adsorbate old_adsorbate = adsorbates[idx_move];

                    // get perturbation of the position in Cartesian space
                    boo::vector<double> dx(3);
                    for (int i_ = 0; i_ < 3; i_++)
                        dx[i_] = parameters.delta * (uniform01(generator) - 0.5);
                    
                    // translate adsorbate (takes care of periodic boundary conditions)
                    adsorbates[idx_move].translate_by_Cartesian_vector(dx, t_matrix, inv_t_matrix, uc_reps);
                    
                    // compute energy of this guest that we propose to move in its new position
                    double E_gf_new = GuestFrameworkEnergy(adsorbates[idx_move], grid_info, energy_grids);
                    double E_gg_new = GuestGuestEnergy(adsorbates, idx_move);
                    double E_new = E_gf_new + E_gg_new;

                    // accept if, loosely, energeticall favorable
                    if (rand_for_acceptance < exp(-(E_new - E_old) / parameters.T)) {
                        stats.N_moves += 1; 
                        E_gg_this_cycle += E_gg_new - E_gg_old; 
                        E_gf_this_cycle += E_gf_new - E_gf_old;
                        if (E_new > 1e6)
                            std::cout << "Move accepted with huge energy" << std::endl;
                        // already overwrote coords with new coords, no need to update coords in this case
                    }
                    else {
                        // replace newly proposed adsorbate position with old position
                        adsorbates[idx_move] = old_adsorbate;
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

                    // pick which adsorbate of this type to change identity of
                    decltype(uniformint.param()) new_range(0, N_g[which_type] - 1); // set new range for rng
                    uniformint.param(new_range);
                    
                    // randomly select guest of this type to propose to change identity of
                    int which_of_this_adsorbate_type = uniformint(generator);  // which guest to propose to delete of this type
                    
                    // loop through adsorbates vector until we find the which_of_this_adsorbate_type'th guest of this type
                    int count_this_type = 0;
                    int idx_id_change;  // index of adsorbate in adsorbates we propose to change identity of 
                    for (idx_id_change = 0; idx_id_change < adsorbates.size(); idx_id_change++) {
                        if (adsorbates[idx_id_change].type == which_type) {
                            if (count_this_type == which_of_this_adsorbate_type)
                                break;
                            count_this_type += 1;
                        }
                    }
                    
                    // compute energy of this guest that we propose to change ID of
                    double E_gf_old = GuestFrameworkEnergy(adsorbates[idx_id_change], grid_info, energy_grids);
                    double E_gg_old = GuestGuestEnergy(adsorbates, idx_id_change);
                    double E_old = E_gf_old + E_gg_old;

                    // randomly pick type to change this adsorbate to. 
                    int change_to_type = which_type;
                    while (change_to_type == which_type)
                        change_to_type = type_generator(generator);
                    
                    if (parameters.debugmode)
                        printf("\tProposal to change identity of guest molecule %d, of type %d, to type %d\n", 
                                idx_id_change, adsorbates[idx_id_change].type, change_to_type);
                    
                    // store old type config so we can replace it later if MC move rejected
                    Adsorbate old_adsorbate = adsorbates[idx_id_change];
                    
                    // insert a new adsorbate at this position
                    // add new guest of this type from the template (easier to start over)
                    Adsorbate new_adsorbate = adsorbatetemplates[change_to_type];
                    // if >1 beads, perform a random rotation
                    if (new_adsorbate.nbeads > 1)
                        PerformUniformRandomRotation(new_adsorbate, generator, std_normal_distn);

                    // translate adsorbate to share same location as first bead
                    boo::vector<double> x_first_bead(3);
                    x_first_bead[0] = adsorbates[idx_id_change].bead_xyz(0, 0);
                    x_first_bead[1] = adsorbates[idx_id_change].bead_xyz(1, 0);
                    x_first_bead[2] = adsorbates[idx_id_change].bead_xyz(2, 0);
                    new_adsorbate.translate_by_Cartesian_vector(x_first_bead, t_matrix, inv_t_matrix, uc_reps);

                    // Finally, replace adsorbate idx_id_change with this new_adsorbate of different identity
                    adsorbates[idx_id_change] = new_adsorbate;
                    
                    // compute energy of this new guest of different identity
                    double E_gf_new = GuestFrameworkEnergy(adsorbates[idx_id_change], grid_info, energy_grids);
                    double E_gg_new = GuestGuestEnergy(adsorbates, idx_id_change);
                    double E_new = E_gf_new + E_gg_new;

                    // Accept move if, loosely, particle identity change was energetically favorable
                    double prob_acceptance_ID_swap = exp(-(E_new - E_old) / parameters.T) * fugacity[change_to_type] / fugacity[which_type] * static_cast<double>(N_g[which_type]) / (static_cast<double>( N_g[change_to_type]) + 1.0);
                    if (rand_for_acceptance < prob_acceptance_ID_swap) { 
                        stats.N_ID_swaps += 1;
                        
                        // keep track of energies
                        E_gg_this_cycle += E_gg_new - E_gg_old;
                        E_gf_this_cycle += E_gf_new - E_gf_old;
                        
                        // update particle numbers
                        N_g[which_type] -= 1; // one less of which type
                        N_g[change_to_type] += 1; // one more of new type
                    }  // end "if we accept this move"
                    else { 
                        // if didn't accept, replace adsorbate with original one
                        adsorbates[idx_id_change] = old_adsorbate;
                    } 
                } // end if N_g > 0
            } // end particle identity swap

            //
            //  MC trial: Regrow (essentially a move)
            //  Delete adsorbate, insert at new position
            //  But in practice, change coords of adsorbate to new position in random location in unit cell
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
                    
                    // randomly select guest of this type to propose to regrow
                    int which_of_this_adsorbate_type = uniformint(generator);  // which guest to propose to delete of this type

                    // loop through adsorbates vector until we find the which_of_this_adsorbate_type'th guest of this type
                    int count_this_type = 0;
                    int idx_regrow;  // index of adsorbate in adsorbates we propose to move
                    for (idx_regrow = 0; idx_regrow < adsorbates.size(); idx_regrow++) {
                        if (adsorbates[idx_regrow].type == which_type) {
                            if (count_this_type == which_of_this_adsorbate_type)
                                break;
                            count_this_type += 1;
                        }
                    }
                      
                    if (parameters.debugmode)
                        printf("\tProposal to regrow guest molecule %d\n", idx_regrow);

                    // compute energy of this guest that we propose to regrow
                    double E_gf_old = GuestFrameworkEnergy(adsorbates[idx_regrow], grid_info, energy_grids);
                    double E_gg_old = GuestGuestEnergy(adsorbates, idx_regrow);
                    double E_old = E_gf_old + E_gg_old;
                    
                    // store this adsorbate's attributes before it was moved
                    // This way, if move is rejected, we can just replace it
                    Adsorbate old_adsorbate = adsorbates[idx_regrow];

                    // insert at new position
                    // add new guest of this type from the template (easier to start over)
                    Adsorbate regrown_adsorbate = adsorbatetemplates[which_type];
                    // if >1 beads, perform a random rotation
                    if (regrown_adsorbate.nbeads > 1)
                        PerformUniformRandomRotation(regrown_adsorbate, generator, std_normal_distn);

                    // translate to these fractional coords
                    boo::vector<double> xf(3);
                    xf[0] =  uniform01(generator) * uc_reps[0];
                    xf[1] =  uniform01(generator) * uc_reps[1];
                    xf[2] =  uniform01(generator) * uc_reps[2];
                    // =  these Cartesian
                    boo::vector<double> x = boo::prod(t_matrix, xf);
                    
                    // translate adsorbate to these coords (also updates fractional coords)
                    regrown_adsorbate.translate_by_Cartesian_vector(x, t_matrix, inv_t_matrix, uc_reps);

                    // replace adsorbate with this one.
                    adsorbates[idx_regrow] = regrown_adsorbate;
                    
                    // compute energy at this new position
                    double E_gf_new = GuestFrameworkEnergy(adsorbates[idx_regrow], grid_info, energy_grids);
                    double E_gg_new = GuestGuestEnergy(adsorbates, idx_regrow);
                    double E_new = E_gf_new + E_gg_new;

                    // accept if, loosely, energeticall favorable
                    if (rand_for_acceptance < exp(-(E_new - E_old) / parameters.T)) {
                        stats.N_regrows += 1; 
                        E_gg_this_cycle += E_gg_new - E_gg_old; E_gf_this_cycle += E_gf_new - E_gf_old;
                        if (E_new > 1e6)
                            std::cout << "Move accepted with huge energy" << std::endl;
                        // already overwrote adsorbate, so no need to do anything
                    }
                    else {
                        // put back old adsorbate
                        adsorbates[idx_regrow] = old_adsorbate;
                    }
                } // end if N_g == 0
            }

            //
            // Collect statistics
            //
            if ((cycle > parameters.numequilibriumtrials) & (cycle_counter % parameters.samplefrequency == 0)) {
                stats.N_samples += 1;
                for (int a_i = 0; a_i < parameters.numadsorbates; a_i++) {
                    stats.N_g_avg[a_i] += N_g[a_i];
                    stats.N_g2_avg[a_i] += N_g[a_i] * N_g[a_i];
                }
                stats.guest_guest_energy_avg += E_gg_this_cycle; // divide by N_samples later
                stats.framework_guest_energy_avg += E_gf_this_cycle; // divide by N_samples later
            }

            //
            // Print stuff if debug mode
            //
            if (parameters.debugmode) {
                // print guest molecule information
                printf("\n\n");
                printf("%lu total guests.\n", adsorbates.size());
                for (int i = 0; i < parameters.numadsorbates; i++) 
                    printf("\tType %d: N_%d = %d\n", i, i, N_g[i]);
                printf("GUEST LIST: (idx, adsorbate type):\n");
                for (int g = 0; g < adsorbates.size(); g++) {
                    printf("(%d, %d) , ", g, adsorbates[g].type);
                }
                printf("\n");
            }
            
            //
            // Make assertions to check code
            //
            if (parameters.makeassertions) {
                double coords_tol = 0.0000001;
                std::vector<int> typecounts(parameters.numadsorbates, 0);
                for (int g = 0; g < adsorbates.size(); g++) {
                    // count types of these beads
                    typecounts[adsorbates[g].type] += 1;

                    // assert first bead is inside unit cell
                    if (OutsideUnitCell(adsorbates[g].bead_xyz_f(0, 0),
                                           adsorbates[g].bead_xyz_f(1, 0),
                                           adsorbates[g].bead_xyz_f(2, 0),
                                             uc_reps) == true) {
                        printf("First bead of adsorbate %d outside unit cell\n", g);
                        printf("uc reps: %d by %d by %d.\n", uc_reps[0], uc_reps[1], uc_reps[2]);
                        for (int bb = 0; bb < adsorbates[g].nbeads; bb++) {
                            printf("Fractional coords bead %d: (%f, %f, %f)\n", bb,
                                                            adsorbates[g].bead_xyz_f(0, bb),
                                                            adsorbates[g].bead_xyz_f(1, bb),
                                                            adsorbates[g].bead_xyz_f(2, bb));
                        }
                        exit(EXIT_FAILURE);
                    }


                    // assert guest beads fractional and cartesian coords match
                    for (int b = 0; b < adsorbates[g].nbeads; b++) {
                        // for each bead
                        boo::matrix_column<boo::matrix<double> > xf_this_bead(adsorbates[g].bead_xyz_f, b);
                        boo::matrix_column<boo::matrix<double> > x_this_bead(adsorbates[g].bead_xyz, b);
                        boo::vector<double> predicted_x = boo::prod(t_matrix, xf_this_bead);
                        if (
                            (predicted_x[0] >  x_this_bead[0] + coords_tol) |
                            (predicted_x[0] <  x_this_bead[0] - coords_tol) |
                            (predicted_x[1] >  x_this_bead[1] + coords_tol) |
                            (predicted_x[1] <  x_this_bead[1] - coords_tol) |
                            (predicted_x[2] >  x_this_bead[2] + coords_tol) |
                            (predicted_x[2] <  x_this_bead[2] - coords_tol)
                        ) {
                            printf("Guest %d, bead %d does not hv match between frac and cart coords.\n", g, b);
                            exit(EXIT_FAILURE);
                        }
                   }  // end loop over beads in guest
                }  // end loop over guests
                int N_g_total = 0;
                for (int i = 0; i < parameters.numadsorbates; i++) {
                    N_g_total += N_g[i];
                    assert(typecounts[i] == N_g[i]);
                }
                assert(N_g_total == adsorbates.size());
                    
            }  // end make assertions
            
//            //
//            // Write adsorbate positions to file (optional)
//            //
//            if (parameters.writeadsorbatepositions) {
//                if ((cycle > parameters.numequilibriumtrials) & (cycle_counter % parameters.writepositionfrequency == 0)) {
//                    N_snapshots ++;
//                    WriteGuestPostionsToFile(adsorbatepositionfile, 
//                              N_g_total,
//                              guestmoleculeinfo,
//                              guestmolecules,
//                              guestbeads,
//                              parameters); 
//                    if (N_snapshots > parameters.num_snapshots) {
//                        printf("Reached %d snapshots, exiting.\n", N_snapshots);
//                        fprintf(outputfile, "\nWrote %d adsorbate snapshot positions in xyz file every %d MC moves.\n", 
//                                        N_snapshots, parameters.writepositionfrequency);
//                        std::exit(EXIT_SUCCESS);
//                    }
//                }
//            }
//
        }  // end inner cycle loop

    }  // end outer cycle loop
    
    // take avg
    stats.guest_guest_energy_avg = 1.0 * stats.guest_guest_energy_avg / stats.N_samples;
    stats.framework_guest_energy_avg = 1.0 * stats.framework_guest_energy_avg / stats.N_samples;

    double sim_time = ReadTimer() - start_of_sim_time;
    
    //
    // check that energy was conserved.
    // Calculate system energy from scratch and compare to energy at this cycle
    //
    fprintf(outputfile, "\nEnergy checks\n");

    double E_gg_system = TotalGuestGuestEnergy(adsorbates);
    double E_gf_system = TotalGuestFrameworkEnergy(adsorbates, grid_info, energy_grids);
    
    fprintf(outputfile, "    E_gg total calc'ed at end: %f K\n", E_gg_system);
    fprintf(outputfile, "    E_gg from adding dE's throughout simulation: %f K\n", E_gg_this_cycle);
    fprintf(outputfile, "    E_gf total calc'ed at end: %f K\n", E_gf_system);
    fprintf(outputfile, "    E_gf from adding dE's throughout simulation: %f K\n", E_gf_this_cycle);
    
    //   
    // write stats to outputfile
    //
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
    
    //
    // write loadings component loadings to outputfile
    //
    for (int n_c = 0; n_c < parameters.numadsorbates; n_c ++) {
        stats.N_g_avg[n_c] = 1.0 * stats.N_g_avg[n_c] / stats.N_samples;
        stats.N_g2_avg[n_c] = 1.0 * stats.N_g2_avg[n_c] / stats.N_samples;
        double N_confidence_bound = sqrt(stats.N_g2_avg[n_c] - stats.N_g_avg[n_c] * stats.N_g_avg[n_c])/sqrt(1.0 * stats.N_samples); // sigma / sqrt(N)

        fprintf(outputfile, "\n    Adsorbate: %s\n", adsorbate[n_c].c_str());
        fprintf(outputfile, "        Fugacity = %f Pa\n", fugacity[n_c]);
        fprintf(outputfile, "        <N_g> (%s) = %f +/- %f molecules/ unit cell\n", adsorbate[n_c].c_str(), 1.0 * stats.N_g_avg[n_c] / uc_reps[0] / uc_reps[1] / uc_reps[2], N_confidence_bound);
        fprintf(outputfile, "        <N_g> (%s) = %f moles/m3\n", adsorbate[n_c].c_str(), stats.N_g_avg[n_c] / volume / 6.022e-7);
        fprintf(outputfile, "        <N_g> (%s) = %f moles/kg framework\n", adsorbate[n_c].c_str(), stats.N_g_avg[n_c] / volume / 6.022e-7 / framework.density);
    }
    fprintf(outputfile, "\n     <E_gg> = %f kJ/mol = %f K\n", stats.guest_guest_energy_avg * 8.314 / 1000.0, stats.guest_guest_energy_avg);
    fprintf(outputfile, "     <E_gf> = %f kJ/mol = %f K", stats.framework_guest_energy_avg * 8.314 / 1000.0, stats.framework_guest_energy_avg);

    fclose(outputfile); 
}
