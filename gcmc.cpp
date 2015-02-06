//
//  GCMC code (mu, V, T) simulations.
//
using namespace std;
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
#include "write_to_outputfile.h"
#include "load_fast_particle_f_array.h"
//#include "pocketblocking.h"
#define min_r .0000000000000001  // don't want to divide by zero...
#define MAX_GUESTS 2000 // max number of guests in an array

void frac_to_cart(double T_matrix[3][3],
    double x_f, double y_f, double z_f,
    double & x, double & y, double & z) 
{ 
    // compute cartesian coordinates from fractional
	x = T_matrix[0][0] * x_f + T_matrix[0][1] * y_f + T_matrix[0][2] * z_f;
	y = T_matrix[1][0] * x_f + T_matrix[1][1] * y_f + T_matrix[1][2] * z_f;
	z = T_matrix[2][0] * x_f + T_matrix[2][1] * y_f + T_matrix[2][2] * z_f;
}

void cart_to_frac(double inv_T_matrix[3][3],
    double & x_f, double & y_f, double & z_f,
    double x, double y, double z)
{ 
    // cartesian to fractional
	x_f = inv_T_matrix[0][0] * x + inv_T_matrix[0][1] * y + inv_T_matrix[0][2] * z;
	y_f = inv_T_matrix[1][0] * x + inv_T_matrix[1][1] * y + inv_T_matrix[1][2] * z;
	z_f = inv_T_matrix[2][0] * x + inv_T_matrix[2][1] * y + inv_T_matrix[2][2] * z;
}

double wrap_0_to_z(double x, double z) {
    // for applying periodic bc's
	return x - z * floor(x / z);
}

//
// Compute potential energy of guest molecule id "which", the contribution from other guests.
//
double E_gg(int N_g, 
    int which, 
    particle_g * guests, 
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
		assert (dx_f < 0.5 * parameters.replication_factor_a);
		assert (dy_f < 0.5 * parameters.replication_factor_b);
		assert (dz_f < 0.5 * parameters.replication_factor_c);
		assert (dx_f > -0.5 * parameters.replication_factor_a);
		assert (dy_f > -0.5 * parameters.replication_factor_b);
		assert (dz_f > -0.5 * parameters.replication_factor_c);
        double dx, dy, dz;
		frac_to_cart(parameters.t_matrix,
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
double E_gg_total(int N_g, 
       particle_g * guests, 
       GCMCParameters parameters)
{
	double E_gg__ = 0.0;
	for (int i = 0; i < N_g; i ++) {
		E_gg__ += E_gg(N_g, i, guests, parameters); 
	}
	return E_gg__ / 2.0; // 2 is for double counting
}

//
//  Get index of flattened, 1D pointer array that corresponds to grid point index (i, j, k)
//
int energy_ind(int i, int j, int k, int N_x, int N_y, int N_z) {
	// returns index in energy_grid corresponding to gridpoint index i,j,k
	return k + j * N_z + i * N_y * N_z;
}

//
// Interpolate energy grid to get guest-framework energy
//  
double E_gf(double x_f_, double y_f_, double z_f_,
            int N_x, int N_y, int N_z,
			double dx_f_, double dy_f_, double dz_f_,
		    double * energy_grid)
{
		// reflect fractional coords in [0,1] for grid
		assert ((x_f_ > 0.0) & (y_f_ > 0.0) & (z_f_ > 0.0));
		assert ((x_f_ < 1.0) & (y_f_ < 1.0) & (z_f_ < 1.0));
		//  FORMAT OF GRID: energy_grid[k+j*N_z+i*N_y*N_z]
		
        // define indices of 8 grid points, lower ones are:
		int i_x_low = floor(x_f_ / dx_f_);
		int i_y_low = floor(y_f_ / dy_f_);
		int i_z_low = floor(z_f_ / dz_f_);
		// trilinear interpolation http://en.wikipedia.org/wiki/Trilinear_interpolation
		
        // difference between our point and the vertices
		double x_d = (x_f_ - i_x_low * dx_f_) / dx_f_;
		double y_d = (y_f_ - i_y_low * dy_f_) / dy_f_;
		double z_d = (z_f_ - i_z_low * dz_f_) / dz_f_;
		
        // smash cube in x direction
		double c00 = energy_grid[energy_ind(i_x_low,i_y_low  ,i_z_low  ,N_x,N_y,N_z)] * (1.0 - x_d) + x_d*energy_grid[energy_ind(i_x_low+1,i_y_low  ,i_z_low  ,N_x,N_y,N_z)];
		double c10 = energy_grid[energy_ind(i_x_low,i_y_low+1,i_z_low  ,N_x,N_y,N_z)] * (1.0 - x_d) + x_d*energy_grid[energy_ind(i_x_low+1,i_y_low+1,i_z_low  ,N_x,N_y,N_z)];
		double c01 = energy_grid[energy_ind(i_x_low,i_y_low  ,i_z_low+1,N_x,N_y,N_z)] * (1.0 - x_d) + x_d*energy_grid[energy_ind(i_x_low+1,i_y_low  ,i_z_low+1,N_x,N_y,N_z)];
		double c11 = energy_grid[energy_ind(i_x_low,i_y_low+1,i_z_low+1,N_x,N_y,N_z)] * (1.0 - x_d) + x_d*energy_grid[energy_ind(i_x_low+1,i_y_low+1,i_z_low+1,N_x,N_y,N_z)];
		
        // further smash cube in y direction
		double c0 = c00 * (1.0 - y_d) + c10 * y_d;
		double c1 = c01 * (1.0 - y_d) + c11 * y_d;
		
        // finally, linear interpolation in z direction
		return c0 * (1 - z_d) + c1 * z_d;
}

//
// Calculate total system framework-guest energy
//
double E_gf_total(int N_g, 
        double * energy_grid0, 
        double * energy_grid1, 
        particle_g * guests, 
        int * N_x, int * N_y, int * N_z, 
        double * dx_f, double * dy_f, double * dz_f)
{
	double E_gf__ = 0.0;
	for (int i = 0; i < N_g; i ++) {
		if (guests[i].type == 0)
			E_gf__ += E_gf(wrap_0_to_z(guests[i].x_f,1.0), wrap_0_to_z(guests[i].y_f,1.0), wrap_0_to_z(guests[i].z_f,1.0),
											N_x[0], N_y[0], N_z[0],
											dx_f[0], dy_f[0], dz_f[0], energy_grid0);
		if (guests[i].type == 1)
			E_gf__ += E_gf(wrap_0_to_z(guests[i].x_f,1.0), wrap_0_to_z(guests[i].y_f,1.0), wrap_0_to_z(guests[i].z_f,1.0),
											N_x[1], N_y[1], N_z[1],
											dx_f[1], dy_f[1], dz_f[1], energy_grid1);
	}
	return E_gf__;
}

//
// Write guest positions to file
//
void write_guest_postions_to_file(particle_g * guests, int N_g, GCMCParameters parameters) {
	FILE * position_file;
	char filename[512] = "output_files/adsorbate_positions.xyz";
	position_file = fopen(filename, "w");
	fprintf(position_file, "%d\n\n", N_g);
	for (int i = 0; i < N_g ; i++) {
		fprintf(position_file, "%s %f %f %f\n", 
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

    readsimulationinputfile(parameters);
    if (parameters.verbose) printf("Read simulation.input\n");
    
    triple_int uc_reps = readunitcellreplicationfile(parameters.frameworkname);
    parameters.replication_factor_a = uc_reps.arg1;
    parameters.replication_factor_b = uc_reps.arg2;
    parameters.replication_factor_c = uc_reps.arg3;
    if (parameters.verbose) printf("Read .uc replication file\n");
    
    parameters.adsorbateMW[0] = get_adsorbate_MW(parameters.adsorbate[0]);
    if (parameters.numadsorbates == 2)
        parameters.adsorbateMW[1] = get_adsorbate_MW(parameters.adsorbate[1]);

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
    
    //
    // Get adsorbate LJ params and store in epsilon and sigma_squared matrices
    //
    pair_double eps_sig = get_guest_FF_params_from_Forcefield(forcefield, parameters.adsorbate[0]);
    parameters.epsilon_matrix[0][0] = eps_sig.arg1; 
    parameters.sigma_squared_matrix[0][0] = eps_sig.arg2 * eps_sig.arg2;
    if (parameters.numadsorbates == 2) {
        eps_sig = get_guest_FF_params_from_Forcefield(forcefield, parameters.adsorbate[1]);
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
	int * N_x = (int *) calloc(2,sizeof(int)); int * N_y = (int *) calloc(2,sizeof(int));	int * N_z = (int *) calloc(2,sizeof(int));
	double * dx_f = (double *) calloc(2,sizeof(double)); double * dy_f = (double *) calloc(2,sizeof(double));	double * dz_f = (double *) calloc(2,sizeof(double));
	double * energy_grid0; double * energy_grid1;
	bool pocket_block_verbose_mode = false;
	for (int n_c = 0; n_c < parameters.numadsorbates; n_c ++) {
		if (parameters.verbose) printf("Importing energy grid %d\n", n_c);
		char gridfilename[512];
        sprintf(gridfilename, "data/grids/%s_%s_%s.txt", framework.name.c_str(), parameters.adsorbate[n_c].c_str(), forcefield.name.c_str());

		ifstream gridfile(gridfilename); // input file stream
		if (gridfile.fail()) {
            printf("Grid file %s could not be loaded...\n", gridfilename);
			exit(EXIT_FAILURE);
		}
		string line;
		getline(gridfile, line);
		istringstream this_line(line);
		this_line >> N_x[n_c] >> N_y[n_c] >> N_z[n_c];
		// create grid import format
		if (n_c == 0)
			energy_grid0 = (double *) calloc(N_x[0] * N_y[0] * N_z[0], sizeof(double));
		if (n_c == 1)
			energy_grid1 = (double *) calloc(N_x[1] * N_y[1] * N_z[1], sizeof(double));

		int i = 0; int j = 0;
		while(getline(gridfile,line)) {
			istringstream this_line(line);
			for (int k = 0; k < N_z[n_c]; k ++) {  // each line is a pencil of z's
				int index_here = k + j * N_z[n_c] + i * N_y[n_c] * N_z[n_c];
				if (n_c == 0)
					this_line >> energy_grid0[index_here];
				if (n_c == 1)
					this_line >> energy_grid1[index_here];
			}
			j += 1; // each line is a pencil of z's for a particular y... so update y index
			if (j == N_y[n_c]) {
				i += 1; // at the end of the z-y sheet, we go to a new x.
				j = 0; // j goes back to zero
			}
		}
		assert(i == N_x[n_c]); // assert that we are at the end of the file
		assert(j == 0);
		gridfile.close();
		// grid spacing in fractional coordinates  (how it was done in writegrid)
		dx_f[n_c] = 1.0/(N_x[n_c] - 1); // grid spacings
		dy_f[n_c] = 1.0/(N_y[n_c] - 1);
		dz_f[n_c] = 1.0/(N_z[n_c] - 1);
		if (parameters.verbose) printf("energy grid %d imported successfully.\n",n_c);
    }
//		//
//		// Flood fill/ pocket blocking
//		//
//		if (parameters.pocket_blocking)
//		{
//			outputfile << "Pocket blocking enabled with temperature " << parameters.T << " K." << endl;
//			double time_before = read_timer();
//			if (pocket_block_verbose_mode == true)
//			{
//				cout << "Pocket blocking beginning. Write a cube for before" << endl;
//				outputfile << "VERBOSE MODE FOR POCKETBLOCKING!" <<endl;
//				if (n_c == 0)
//					write_cube("before_blocking_0", framework, &parameters, energy_grid0, N_x[n_c],N_y[n_c],N_z[n_c]); //just so you can see the grid before ...
//				if (n_c == 1)
//					write_cube("before_blocking_1", framework, &parameters, energy_grid1, N_x[n_c],N_y[n_c],N_z[n_c]); //just so you can see the grid before ...
//				cout << "Cube written." << endl;
//			}
//			int num_pockets = -1;
//			if (n_c == 0)
//				num_pockets = find_and_block_pockets(energy_grid0, N_x[0], N_y[0], N_z[0], parameters.T, parameters);
//			if (n_c == 1)
//				num_pockets = find_and_block_pockets(energy_grid1, N_x[1], N_y[1], N_z[1], parameters.T, parameters);
//			outputfile << parameters.adsorbate[n_c] << " : found " << num_pockets << " inaccessible pockets" << endl;
//			// energy_grid is passed as a poitner, so it is modified now if accessible pockets were there.
//			if (pocket_block_verbose_mode == true)
//			{
//				cout << "Pocket blocking finished. Write a cube for after" << endl;
//				double time_after = read_timer();
//				outputfile << "Time spent to find and block pockets: " << time_after-time_before << endl;
//				if (n_c == 0)
//					write_cube("after_blocking_0", framework, &parameters, energy_grid0, N_x[n_c],N_y[n_c],N_z[n_c]); // ... and after pocket blocking
//				if (n_c == 1)
//					write_cube("after_blocking_1", framework, &parameters, energy_grid1, N_x[n_c],N_y[n_c],N_z[n_c]); // ... and after pocket blocking
//				cout << "Cube written." << endl;
//			}
//		}
//	}
	
    //
    // Initialize stats and guests array (includes both components via "type" attribute)
    //
	GCMC_stats stats;
    initialize_GCMC_stats(stats);
	
    double volume = framework.volume_unitcell * parameters.replication_factor_a * parameters.replication_factor_b * parameters.replication_factor_c; // A ^ 3
	
    particle_g * guests = (particle_g *) malloc(MAX_GUESTS * sizeof(particle_g));
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

	//
	//  RUN GCMC
	//
	double start_of_sim_time=read_timer();
	if (parameters.verbose) printf("Starting simulation...\n");
//
//
//	// initialize
//	double E_gg_this_cycle = 0.0; double E_gf_this_cycle = 0.0; // assumes we start with zero particles
//	int cycle_counter = 0;
//	int adsorbate_index_list[2][MAX_GUESTS] = { -1 }; // keep indices of particles here
//	// for adaptive algorithm for step size
//	int block_no = 0; int N_accepted_translations = 0; int N_attempted_translations = 0;
//	double ratio_old = 0.0; double delta_old = parameters.delta;
//	for (int cycle = 0; cycle < parameters.numtrials; cycle++)
//  {
//		int ninnercycles = (N_g[0]+N_g[1]) > 20 ? (N_g[0]+N_g[1]) : 20; // makes proportional to number of guests
//		for (int inner_cycle = 0 ; inner_cycle < ninnercycles ; inner_cycle ++)
//		{
//			assert(N_g_total == N_g[0] + N_g[1]);
//			cycle_counter += 1;
//			double which_move = uniform01(generator); // particle swap or translation?
//			double rand_for_acceptance = uniform01(generator); // for testing acceptance
//			int which_type = type_generator(generator); // select adsorbate type . e.g. Xe or Kr in binary mixture.
//	//		// testing energy
//	//		N_g = 5;
//	//		double x_f,y_f,z_f;
//	//		guests[0].x = 31.688; guests[0].y =   19.090; guests[0].z =   8.025;
//	//		guests[1].x = 19.790; guests[1].y =   17.562; guests[1].z =  17.791;
//	//		guests[2].x = 14.261; guests[2].y =    6.137; guests[2].z =   5.837;
//	//		guests[3].x = 17.466; guests[3].y =   19.442; guests[3].z =   3.450;
//	//		guests[4].x = 24.163; guests[4].y =   22.185; guests[4].z =  23.055;
//	//		for (int i=0; i < N_g; i++)
//	//		{
//	//			cart_to_frac(parameters.inv_T_matrix,x_f,y_f,z_f,guests[i].x,guests[i].y,guests[i].z);
//	//			guests[i].x_f = x_f;
//	//			guests[i].y_f = y_f;
//	//			guests[i].z_f = z_f;
//	//		}
//	//		double E_gg_total = 0.0;
//	//		for (int i=0; i < N_g; i++)
//	//		{
//	//			E_gg_total += E_gg(N_g, i, guests, parameters)/2.0;
//	//		}
//	//		cout << "E_gg total = " << E_gg_total;
//	//		exit(EXIT_FAILURE);
//			//
//			//  MC TRIAL: INSERTION
//			//
//			if (which_move < parameters.p_exchange/2.0)
//			{
//				stats.N_insertion_trials += 1;
//				// insertion location @ these fractional coordinates
//				guests[N_g_total].x_f = uniform01(generator) * parameters.replication_factor_a;
//				guests[N_g_total].y_f = uniform01(generator) * parameters.replication_factor_b;
//				guests[N_g_total].z_f = uniform01(generator) * parameters.replication_factor_c;
//				double x,y,z; frac_to_cart(parameters.T_matrix, guests[N_g_total].x_f, guests[N_g_total].y_f, guests[N_g_total].z_f, x, y, z);
//				guests[N_g_total].x = x;
//				guests[N_g_total].y = y;
//				guests[N_g_total].z = z;
//				guests[N_g_total].type = which_type; // declare particle type
//				// compute energies here
//				double energy_gf = 1e6;
//				if (which_type == 0)
//					energy_gf = E_gf(wrap_0_to_z(guests[N_g_total].x_f , 1.0),wrap_0_to_z(guests[N_g_total].y_f , 1.0), wrap_0_to_z(guests[N_g_total].z_f ,1.0),
//																											N_x[which_type], N_y[which_type], N_z[which_type],
//																											dx_f[which_type], dy_f[which_type], dz_f[which_type],
//																										energy_grid0);
//				if (which_type == 1)
//					energy_gf = E_gf(wrap_0_to_z(guests[N_g_total].x_f , 1.0),wrap_0_to_z(guests[N_g_total].y_f , 1.0), wrap_0_to_z(guests[N_g_total].z_f ,1.0),
//																											N_x[which_type], N_y[which_type], N_z[which_type],
//																											dx_f[which_type], dy_f[which_type], dz_f[which_type],
//																											energy_grid1);
//				double energy_gg = E_gg(N_g_total, N_g_total, guests, parameters, adsorbate_sigma_matrix, adsorbate_epsilon_matrix);
//				double E_insertion = energy_gf + energy_gg;
//				// accept if energetically favorable
//				double acceptance_insertion = parameters.fugacity[which_type] * volume / ( (N_g[which_type] + 1) * 1.3806488e7 * parameters.T ) * exp(-E_insertion/parameters.T);
//				if (parameters.verbose)
//				{
//					cout << "Adsorbate id list: ";
//					for (int ad = 0 ; ad < N_g[0] ; ad++)
//					 cout << adsorbate_index_list[0][ad] << endl;
//					printf("\tProposal to insert at xf=%f,yf=%f,zf=%f\n.",guests[N_g_total].x_f,guests[N_g_total].y_f,guests[N_g_total].z_f);
//					cout << "\tE_gg = " << energy_gg << endl;
//					cout << "\tE_gf = " << energy_gf << endl;
//					cout << "\tAcceptance prob = " << acceptance_insertion << endl;
//				}
//				if (rand_for_acceptance < acceptance_insertion) // accept insertion move
//				{
//					if (parameters.verbose) cout << "\tInsertion accepted. " << endl;
//					stats.N_insertions += 1;
//						if (energy_gf > 1e6)
//							cout << "Insertion accepted with huge energy" << endl;
//					// add adsorbate guests index to adsorbate_index_list
//					adsorbate_index_list[which_type][N_g[which_type]] = N_g_total;
//					N_g[which_type] += 1;
//					N_g_total += 1;
//					E_gg_this_cycle += energy_gg; E_gf_this_cycle += energy_gf;
//				}
//			}
//			//
//			//  MC TRIAL: DELETION
//			//
//			else if (which_move < parameters.p_exchange)
//			{
//				if (parameters.verbose) cout << "Deletion Trial." << endl;
//				stats.N_deletion_trials += 1;
//				if (N_g[which_type] > 0)
//				{
//					// set new range for uniform int generator.
//					decltype(uniformint.param()) new_range (0, N_g[which_type] - 1);
//					uniformint.param(new_range);
//					int idx_delete_type = uniformint(generator); //mt_lrand() % N_g; // select idx of guest to delete in adsorbate_index_list
//					assert(idx_delete_type < N_g[which_type]);
//					int idx_delete = adsorbate_index_list[which_type][idx_delete_type]; // corresponding global ID in guests
//					if (idx_delete >= N_g_total)
//					{
//						printf("Idx delete %d, N_g %d, Ng total %d,idx delete type %d\n",idx_delete,N_g[which_type],N_g_total,idx_delete_type);
//						printf("adsorbate index list: ");
//						for (int ad = 0 ; ad < N_g[0] ; ad++)
//						 cout << adsorbate_index_list[0][ad] << endl;
//					}
//					assert(idx_delete < N_g_total);
//					// compute energy of these guest
//					double energy_gf = 1e6;
//					if (which_type == 0)
//						energy_gf = E_gf(wrap_0_to_z(guests[idx_delete].x_f , 1.0),wrap_0_to_z(guests[idx_delete].y_f , 1.0), wrap_0_to_z(guests[idx_delete].z_f ,1.0),
//																												N_x[which_type], N_y[which_type], N_z[which_type],
//																												dx_f[which_type], dy_f[which_type], dz_f[which_type],
//																												energy_grid0);
//					if (which_type == 1)
//						energy_gf = E_gf(wrap_0_to_z(guests[idx_delete].x_f , 1.0),wrap_0_to_z(guests[idx_delete].y_f , 1.0), wrap_0_to_z(guests[idx_delete].z_f ,1.0),
//																												N_x[which_type], N_y[which_type], N_z[which_type],
//																												dx_f[which_type], dy_f[which_type], dz_f[which_type],
//																												energy_grid1);
//					double energy_gg = E_gg(N_g_total, idx_delete, guests, parameters, adsorbate_sigma_matrix, adsorbate_epsilon_matrix);
//					double E_deletion = energy_gf + energy_gg;
//					double acceptance_del = ( N_g[which_type] * 1.3806488e7 * parameters.T ) / (parameters.fugacity[which_type] * volume) * exp(E_deletion/parameters.T);
//					if (rand_for_acceptance < acceptance_del)
//					{
//						if (energy_gf > 1e6)
//						{
//							cout << "Deletion accepted with huge energy" << endl;
//							printf("N_g = %d, energy_gg = %f, idx_delete = %d, idx_delete_type = %d, N_g1 = %d\n", N_g[0],energy_gg,idx_delete, idx_delete_type,N_g[1]);
//						}
//						stats.N_deletions += 1;
//						// move last particle to here
//						guests[idx_delete] = guests[N_g_total - 1];
////						guests[idx_delete].x_f = guests[N_g_total - 1].x_f;	guests[idx_delete].y_f = guests[N_g_total - 1].y_f; guests[idx_delete].z_f = guests[N_g_total - 1].z_f;
////						guests[idx_delete].x = guests[N_g_total - 1].x;	guests[idx_delete].y = guests[N_g_total - 1].y; guests[idx_delete].z = guests[N_g_total - 1].z;
////						guests[idx_delete].type = guests[N_g_total-1].type;
//						// if delete Xe in the middle of adsorbate_index_list, move the last one here so we see it.
//						adsorbate_index_list[which_type][idx_delete_type] = adsorbate_index_list[which_type][N_g[which_type] - 1];
//						// also we moved the last guest to the middle... find this!
//						for (int orz = 0 ; orz < N_g[guests[N_g_total - 1].type] ; orz ++)
//						{
//							if (adsorbate_index_list[guests[N_g_total-1].type][orz] == N_g_total - 1)
//							{
//								adsorbate_index_list[guests[N_g_total-1].type][orz] = idx_delete; //now this one is at idx_delete
//								break;
//							}
//						}
//						N_g[which_type] -= 1;
//						N_g_total -= 1;
//						E_gg_this_cycle -= energy_gg; E_gf_this_cycle -= energy_gf;
//						if (parameters.verbose)
//						{
//							cout << "Deletion accepted with probability " << acceptance_del << endl;
//						}
//					}
//				}
//			}
//			//
//			//	MC TRIAL:  MOVE
//			//
//			else if (which_move < parameters.p_exchange + parameters.p_move)
//			{
//				if (parameters.verbose) cout << "Translation Trial." << endl;
//				stats.N_move_trials += 1;
//				if (N_g[which_type] > 0)
//				{
//					N_attempted_translations += 1;
//					// Randomly choose an adsorabte of which_type
//					decltype(uniformint.param()) new_range (0, N_g[which_type]-1); // set new range for rng
//					uniformint.param(new_range);
//					int idx_move_type = uniformint(generator);
//					assert(idx_move_type < N_g[which_type]);
//					int idx_move = adsorbate_index_list[which_type][idx_move_type]; // global ID in guests
//					assert(idx_move < N_g_total);
//		//			if (idx_move == N_g - 1) cout << "found last one" << endl;
//		//			if (idx_move == 0) cout << "found first one" << endl;
//					// get old energy
//					double energy_gf_old = 0.0;
//					if (which_type == 0)
//						energy_gf_old = E_gf(wrap_0_to_z(guests[idx_move].x_f , 1.0),wrap_0_to_z(guests[idx_move].y_f , 1.0), wrap_0_to_z(guests[idx_move].z_f ,1.0),
//																												N_x[which_type], N_y[which_type], N_z[which_type],
//																												dx_f[which_type], dy_f[which_type], dz_f[which_type],
//																												energy_grid0);
//					if (which_type == 1)
//						energy_gf_old = E_gf(wrap_0_to_z(guests[idx_move].x_f , 1.0),wrap_0_to_z(guests[idx_move].y_f , 1.0), wrap_0_to_z(guests[idx_move].z_f ,1.0),
//																												N_x[which_type], N_y[which_type], N_z[which_type],
//																												dx_f[which_type], dy_f[which_type], dz_f[which_type],
//																												energy_grid1);
//					double energy_gg_old = E_gg(N_g_total, idx_move, guests, parameters, adsorbate_sigma_matrix, adsorbate_epsilon_matrix);
//					// store old coords
//					particle_g old_position = guests[idx_move];
//					// perturb by this amount.
//					guests[idx_move].x += parameters.delta * (uniform01(generator) - 0.5);
//					guests[idx_move].y += parameters.delta * (uniform01(generator) - 0.5);
//					guests[idx_move].z += parameters.delta * (uniform01(generator) - 0.5);
//					// convert back to frac
//					double x_f, y_f, z_f;
//					cart_to_frac(parameters.inv_T_matrix,x_f,y_f,z_f,guests[idx_move].x,guests[idx_move].y,guests[idx_move].z);
//					guests[idx_move].x_f = x_f;
//					guests[idx_move].y_f = y_f;
//					guests[idx_move].z_f = z_f;
//					// if move outside of box...
//					if ((x_f > parameters.replication_factor_a) | (y_f > parameters.replication_factor_b) | (z_f > parameters.replication_factor_c) | (x_f < 0.0) | (y_f < 0.0) | (z_f < 0.0) )
//					{
//						// adjust fractional coordinates, then later Cartesian coords
//						guests[idx_move].x_f = wrap_0_to_z(x_f,parameters.replication_factor_a);
//						guests[idx_move].y_f = wrap_0_to_z(y_f,parameters.replication_factor_b);
//						guests[idx_move].z_f = wrap_0_to_z(z_f,parameters.replication_factor_c);
//						double x,y,z;
//						frac_to_cart(parameters.T_matrix,guests[idx_move].x_f,guests[idx_move].y_f,guests[idx_move].z_f,x,y,z);
//						guests[idx_move].x = x;
//						guests[idx_move].y = y;
//						guests[idx_move].z = z;
//					}
//					assert(guests[idx_move].x_f > 0.0); assert(guests[idx_move].x_f < parameters.replication_factor_a);
//					assert(guests[idx_move].y_f > 0.0); assert(guests[idx_move].y_f < parameters.replication_factor_b);
//					assert(guests[idx_move].z_f > 0.0); assert(guests[idx_move].z_f < parameters.replication_factor_c);
//					// get new energy
//					double energy_gf_new = 1e6;
//					if (which_type == 0)
//						energy_gf_new = E_gf(wrap_0_to_z(guests[idx_move].x_f , 1.0),wrap_0_to_z(guests[idx_move].y_f , 1.0), wrap_0_to_z(guests[idx_move].z_f ,1.0),
//																												N_x[which_type], N_y[which_type], N_z[which_type],
//																												dx_f[which_type], dy_f[which_type], dz_f[which_type],
//																												energy_grid0);
//					if (which_type == 1)
//						energy_gf_new = E_gf(wrap_0_to_z(guests[idx_move].x_f , 1.0),wrap_0_to_z(guests[idx_move].y_f , 1.0), wrap_0_to_z(guests[idx_move].z_f ,1.0),
//																												N_x[which_type], N_y[which_type], N_z[which_type],
//																												dx_f[which_type], dy_f[which_type], dz_f[which_type],
//																												energy_grid1);
//					double energy_gg_new = E_gg(N_g_total, idx_move, guests, parameters, adsorbate_sigma_matrix, adsorbate_epsilon_matrix);
//
//					double dE = energy_gf_new + energy_gg_new - (energy_gf_old + energy_gg_old);
//					if (rand_for_acceptance < exp(-dE/parameters.T))
//					{
//						stats.N_moves += 1; N_accepted_translations += 1;
//						E_gg_this_cycle += energy_gg_new - energy_gg_old; E_gf_this_cycle += energy_gf_new - energy_gf_old;
//						if (energy_gf_new > 1e6)
//						{
//							cout << "Move accepted with huge energy" << endl;
//							cout << "Old energy = " << energy_gf_old << endl;
//							cout << "New energy = " << energy_gf_new << endl;
//						}
//						// already overwrote coords with new coords
//					}
//					else
//					{
//						// replace new cords with old coords
//						guests[idx_move] = old_position;
//					}
//				} // end if N_g == 0
//			} // end translation
//
//			//
//			//  PARTICLE IDENTITY CHANGE FOR DUAL COMPONENT
//			//
//			else
//			{
//				if (N_g[which_type] > 0) // if there are paricles of this type in the system
//				{
//					stats.N_ID_swap_trials += 1;
//					// pick which particles to change identity of
//					decltype(uniformint.param()) new_range (0, N_g[which_type] - 1); // set new range for rng
//					uniformint.param(new_range);
//					int idx_type = uniformint(generator); // which of this component?
//					int idx = adsorbate_index_list[which_type][idx_type]; // global index of this component
//					assert(guests[idx].type == which_type);
//
//					// get energy difference with current ID
//					double energy_gg_old = E_gg(N_g_total, idx , guests, parameters, adsorbate_sigma_matrix, adsorbate_epsilon_matrix);
//					double energy_gf_0 = E_gf(wrap_0_to_z(guests[idx].x_f , 1.0),wrap_0_to_z(guests[idx].y_f , 1.0), wrap_0_to_z(guests[idx].z_f ,1.0),
//																											N_x[0], N_y[0], N_z[0],
//																											dx_f[0], dy_f[0], dz_f[0],
//																										energy_grid0);
//					double energy_gf_1 = E_gf(wrap_0_to_z(guests[idx].x_f , 1.0),wrap_0_to_z(guests[idx].y_f , 1.0), wrap_0_to_z(guests[idx].z_f ,1.0),
//																											N_x[1], N_y[1], N_z[1],
//																											dx_f[1], dy_f[1], dz_f[1],
//																										energy_grid1);
//					// change identity, recalculate energy
//					if (which_type == 0)
//						guests[idx].type = 1;
//					if (which_type == 1)
//						guests[idx].type = 0;
//					double energy_gg_new = E_gg(N_g_total, idx , guests, parameters, adsorbate_sigma_matrix, adsorbate_epsilon_matrix);
//					double dE;
//					if (which_type == 0) // if changing it to a 1
//						dE = energy_gf_1 + energy_gg_new - (energy_gf_0 + energy_gg_old);
//					if (which_type == 1) // if changing it to a 0
//						dE = energy_gf_0 + energy_gg_new - (energy_gf_1 + energy_gg_old);
//					int new_type = guests[idx].type;
//					double prob_acceptance_ID_swap = exp( - dE/ parameters.T) * parameters.fugacity[new_type] / parameters.fugacity[which_type] * (double) N_g[which_type] / ((double) N_g[new_type]+ 1.0);
//					if (rand_for_acceptance < prob_acceptance_ID_swap) // ACCEPT
//					{
//						stats.N_ID_swaps += 1;
//						// keep guest particle as new declared type.
//						// keep track of energies
//						E_gg_this_cycle += energy_gg_new - energy_gg_old;
//						if (which_type == 0) // if changed to a 1
//							E_gf_this_cycle += energy_gf_1 - energy_gf_0;
//						if (which_type == 1) // if changed to a zero
//							E_gf_this_cycle += energy_gf_0 - energy_gf_1;
//						// deal with adorbate ID list
//						adsorbate_index_list[which_type][idx_type] = adsorbate_index_list[which_type][N_g[which_type] - 1]; // move last one of which_type to here, this one is deleted.
//						adsorbate_index_list[new_type][N_g[new_type]] = idx; // now this index is part of the other component ID
//						// update particle numbers
//						N_g[which_type] -= 1; // one less of which type
//						N_g[new_type] += 1; // one more of new type
//					}
//					else // if didn't accept
//						guests[idx].type = which_type; // go back to old type.
//				} // end if N_g > 0
//			} // end particle identity swap
//
//			// assert N_g < MAX_GUESTS
//			if (N_g_total >= MAX_GUESTS - 2)
//			{
//				cerr << "N_g > MAX_GUESTS!!! " << endl;
//				outputfile << "N_g > MAX_GUESTS!!! " << endl;
//				exit(EXIT_FAILURE);
//			}
//
//			//
//			// ADAPT delta to try to get 20% acceptance moves
//			//
////			if ((cycle < parameters.equilibrium_trials) & (N_attempted_translations % 1000 == 0) & (N_attempted_translations != 0))
////			{
////				cout << "N attempt " << N_attempted_translations << endl;
////				cout << "N accept " << N_accepted_translations << endl;
////				cout << "Acceptance ratio = " << 1.0 * N_accepted_translations / N_attempted_translations << endl;
////				cout << "Delta before = " << parameters.delta << endl;
//////				adapt_delta(ratio_old, delta_old, N_attempted_translations, N_accepted_translations, parameters, block_no);
////				block_no += 1;
////				cout << "Delta after = " << parameters.delta << endl;
////			}
//			//
//			// COLLECT STATS
//			//
//			if ((cycle > parameters.equilibrium_trials) & (cycle_counter % parameters.sample_every == 0))
//			{
//				stats.N_samples += 1;
//				stats.N_g_avg[0] += N_g[0]; stats.N_g_avg[1] += N_g[1];
//				stats.N_g2_avg[0] += N_g[0] * N_g[0]; stats.N_g2_avg[1] += N_g[1] * N_g[1];
//				stats.guest_guest_energy_avg += E_gg_this_cycle;
//				stats.framework_guest_energy_avg += E_gf_this_cycle;
//			}
//			if (parameters.verbose == true)
//			{
//				cout << "MC move " << cycle_counter << endl;
//				cout << "\tN_g [0]= " << N_g[0] << endl;
//				cout << "\tN_g [1]= " << N_g[1] << endl;
//				cout << "\tN_g_total = " << N_g_total << endl;
//				cout << "\tE_gg = " << E_gg_this_cycle << endl;
//				cout << "\tE_gf = " << E_gf_this_cycle << endl;
//			}
//		} // end inner cycle
//
//	}
//
//	stats.guest_guest_energy_avg = 1.0 * stats.guest_guest_energy_avg / stats.N_samples;
//	stats.framework_guest_energy_avg = 1.0 * stats.framework_guest_energy_avg / stats.N_samples;
//
//	// write positions of adosrbates to file
////	write_guest_postions_to_file(guests,N_g[0]+N_g[1],parameters);
//
// 	double sim_time=read_timer()-start_of_sim_time;
//	outputfile << " sim time (s): " << sim_time<<endl;
//	outputfile << "-----Energy checks.------" << endl;
//	// check energy calcs
//	double E_gg_system = E_gg_total(N_g[0]+N_g[1], guests, parameters, adsorbate_sigma_matrix, adsorbate_epsilon_matrix);
//	double E_gf_system = E_gf_total(N_g[0]+N_g[1], energy_grid0, energy_grid1, guests, N_x,N_y,N_z, dx_f,dy_f,dz_f);
//	outputfile << "\tE_gg total = " << E_gg_system << endl;
//	outputfile << "\tE_gg total, using MC = " << E_gg_this_cycle << endl;
//	outputfile << "\tE_gf total = " << E_gf_system << endl;
//	outputfile << "\tE_gf total, using MC = " << E_gf_this_cycle << endl;
//
//	outputfile << "\n\nSTATS\n=======\n";
//	outputfile << "N_insertions: " << stats.N_insertions << " / " << stats.N_insertion_trials << endl;
//	outputfile << "N_deletions: " << stats.N_deletions << " / " << stats.N_deletion_trials << endl;
//	outputfile << "N_moves: " << stats.N_moves << " / " << stats.N_move_trials << endl;
//	outputfile << "N_ID_Swaps: " << stats.N_ID_swaps << " / " << stats.N_ID_swap_trials << endl;
//	outputfile << "N_samples: " << stats.N_samples << endl;
//
//	outputfile << "\n\nRESULTS\n=======\n";
//	outputfile << "\tTemperature (K): " << parameters.T << endl;
//	for (int n_c = 0 ; n_c < parameters.numadsorbates ; n_c ++)
//	{
//		string adsorbate;
//		if (n_c == 0)
//			adsorbate = parameters.adsorbate[0];
//		if (n_c == 1)
//			adsorbate = parameters.adsorbate[1];
//		stats.N_g_avg[n_c] = 1.0 * stats.N_g_avg[n_c] / stats.N_samples;
//		stats.N_g2_avg[n_c] = 1.0 * stats.N_g2_avg[n_c] / stats.N_samples;
//		double N_confidence_bound = sqrt(stats.N_g2_avg[n_c] - stats.N_g_avg[n_c] * stats.N_g_avg[n_c])/sqrt(1.0 * stats.N_samples); // sigma / sqrt(N)
//		outputfile << "Adsorbate " << n_c << ": " << adsorbate <<endl;
//		outputfile << "\tFugacity (Pa) = " << parameters.fugacity[n_c] << endl;
//		outputfile << "\t<N_g> (" << adsorbate << ") = " << 1.0 * stats.N_g_avg[n_c] / parameters.replication_factor_a / parameters.replication_factor_b / parameters.replication_factor_c << " +/- " << N_confidence_bound << " molecules per unit cell" << endl;
//		outputfile << "\t<N_g> (" << adsorbate << ", moles/m3) = " << stats.N_g_avg[n_c]/volume/6.022e-7 << endl;
//		outputfile << "\t<N_g> (" << adsorbate << ", moles/kg framework) = " << stats.N_g_avg[n_c]/volume/6.022e-7/framework.density << endl;
//	}
//	outputfile << "<E_gg> (K) = " << stats.guest_guest_energy_avg <<endl; //* 8.314 / 1000.0<< endl;
//	outputfile << "\t(kJ/mol) = " << stats.guest_guest_energy_avg * 8.314 / 1000.0<< endl;
//	outputfile << "<E_gf> (K) = " << stats.framework_guest_energy_avg <<endl; //* 8.314 / 1000.0<< endl;
//	outputfile << "\t(kJ/mol) = " << stats.framework_guest_energy_avg * 8.314 / 1000.0<< endl;
}
