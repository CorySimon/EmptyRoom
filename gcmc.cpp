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
//double E_gg(int N_g, 
//    int which, 
//    particle_g * guests, 
//    Parameters parameters, 
//    double adsorbate_sigma_matrix[2][2], 
//    double adsorbate_epsilon_matrix[2][2])
//{
//	double E_gg_ = 0.0; // initiate
//	for (int k = 0 ; k < N_g; k++ ) { 
//        // nearest particle image
//		if (k == which) 
//            continue; // do not include self interation orz
//		// use nearest image
//		double dx_f = guests[k].x_f - guests[which].x_f;
//		double dy_f = guests[k].y_f - guests[which].y_f;
//		double dz_f = guests[k].z_f - guests[which].z_f;
//		dx_f = dx_f - parameters.replication_factor_a * round(dx_f / parameters.replication_factor_a);
//		dy_f = dy_f - parameters.replication_factor_b * round(dy_f / parameters.replication_factor_b);
//		dz_f = dz_f - parameters.replication_factor_c * round(dz_f / parameters.replication_factor_c);
//		assert (dx_f < 0.5 * parameters.replication_factor_a);
//		assert (dy_f < 0.5 * parameters.replication_factor_b);
//		assert (dz_f < 0.5 * parameters.replication_factor_c);
//		assert (dx_f > -0.5 * parameters.replication_factor_a);
//		assert (dy_f > -0.5 * parameters.replication_factor_b);
//		assert (dz_f > -0.5 * parameters.replication_factor_c);
//        double dx, dy, dz;
//		frac_to_cart(parameters.T_matrix,
//            dx_f, dy_f, dz_f,
//            dx, dy, dz);
//		double r2 = dx*dx + dy*dy + dz*dz;
//		if (r2 < min_r) 
//            return 100000000000000.0; // overwrite if too small
//		if (r2 < parameters.r_cutoff_squared) {
//			double sigma_over_r_sixth = pow(adsorbate_sigma_matrix[guests[which].type][guests[k].type] * adsorbate_sigma_matrix[guests[which].type][guests[k].type] /r2,3.0); // TODO: make sig2
//			E_gg_ += 4.0 * adsorbate_epsilon_matrix[guests[which].type][guests[k].type] * sigma_over_r_sixth * (sigma_over_r_sixth - 1.0); // energy of framework atom id with guest_atom
//		}
//	}
//	return E_gg_;
//}

//
// calculate total system guest-guest energy
// 
//double E_gg_total(int N_g, 
//        particle_g * guests, 
//        Parameters parameters, 
//        double adsorbate_sigma_matrix[2][2], 
//        double adsorbate_epsilon_matrix[2][2])
//{
//	double E_gg__ = 0.0;
//	for (int i = 0; i < N_g ; i ++) {
//		E_gg__ += E_gg(N_g, i, guests, parameters, adsorbate_sigma_matrix, adsorbate_epsilon_matrix)/2.0; // 2 is for double counting
//	}
//	return E_gg__;
//}

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
//double E_gf(double x_f_, double y_f_, double z_f_,
//            int N_x, int N_y, int N_z,
//			double dx_f_, double dy_f_, double dz_f_,
//		    double * energy_grid)
//
//{
//		// reflect fractional coords in [0,1] for grid
//		assert ((x_f_ > 0.0) & (y_f_ > 0.0) & (z_f_ > 0.0));
//		assert ((x_f_ < 1.0) & (y_f_ < 1.0) & (z_f_ < 1.0));
//		//  FORMAT OF GRID: energy_grid[k+j*N_z+i*N_y*N_z]
//		
//        // define indices of 8 grid points, lower ones are:
//		int i_x_low = floor(x_f_/dx_f_);
//		int i_y_low = floor(y_f_/dy_f_);
//		int i_z_low = floor(z_f_/dz_f_);
//		// trilinear interpolation http://en.wikipedia.org/wiki/Trilinear_interpolation
//		
//        // difference between our point and the vertices
//		double x_d = (x_f_-i_x_low*dx_f_)/dx_f_;
//		double y_d = (y_f_-i_y_low*dy_f_)/dy_f_;
//		double z_d = (z_f_-i_z_low*dz_f_)/dz_f_;
//		
//        // smash cube in x direction
//		double c00 = energy_grid[energy_ind(i_x_low,i_y_low  ,i_z_low  ,N_x,N_y,N_z)] * (1.0 - x_d) + x_d*energy_grid[energy_ind(i_x_low+1,i_y_low  ,i_z_low  ,N_x,N_y,N_z)];
//		double c10 = energy_grid[energy_ind(i_x_low,i_y_low+1,i_z_low  ,N_x,N_y,N_z)] * (1.0 - x_d) + x_d*energy_grid[energy_ind(i_x_low+1,i_y_low+1,i_z_low  ,N_x,N_y,N_z)];
//		double c01 = energy_grid[energy_ind(i_x_low,i_y_low  ,i_z_low+1,N_x,N_y,N_z)] * (1.0 - x_d) + x_d*energy_grid[energy_ind(i_x_low+1,i_y_low  ,i_z_low+1,N_x,N_y,N_z)];
//		double c11 = energy_grid[energy_ind(i_x_low,i_y_low+1,i_z_low+1,N_x,N_y,N_z)] * (1.0 - x_d) + x_d*energy_grid[energy_ind(i_x_low+1,i_y_low+1,i_z_low+1,N_x,N_y,N_z)];
//		
//        // further smash cube in y direction
//		double c0 = c00 * (1.0 - y_d) + c10 * y_d;
//		double c1 = c01 * (1.0 - y_d) + c11 * y_d;
//		
//        // finally, linear interpolation in z direction
//		return c0 * (1 - z_d) + c1 * z_d;
//}
//
////
//// Calculate total system framework-guest energy
////
//double E_gf_total(int N_g, 
//        double * energy_grid0, 
//        double * energy_grid1, 
//        particle_g * guests, 
//        int * N_x, int * N_y, int * N_z, 
//        double * dx_f, double * dy_f, double * dz_f)
//{
//	double E_gf__ = 0.0;
//	for (int i = 0; i < N_g ; i ++) {
//		if (guests[i].type == 0)
//			E_gf__ += E_gf(wrap_0_to_z(guests[i].x_f,1.0), wrap_0_to_z(guests[i].y_f,1.0), wrap_0_to_z(guests[i].z_f,1.0),
//											N_x[guests[i].type] ,N_y[guests[i].type],N_z[guests[i].type],
//											dx_f[guests[i].type],dy_f[guests[i].type],dz_f[guests[i].type], energy_grid0);
//		if (guests[i].type == 1)
//			E_gf__ += E_gf(wrap_0_to_z(guests[i].x_f,1.0), wrap_0_to_z(guests[i].y_f,1.0), wrap_0_to_z(guests[i].z_f,1.0),
//											N_x[guests[i].type] ,N_y[guests[i].type],N_z[guests[i].type],
//											dx_f[guests[i].type],dy_f[guests[i].type],dz_f[guests[i].type], energy_grid1);
//	}
//	return E_gf__;
//}
//
////
//// Write guest positions to file
////
//void write_guest_postions_to_file(particle_g * guests, int N_g, Parameters parameters) {
//	FILE * position_file;
//	char filename[512] = "output_files/adsorbate_positions.xyz";
//	position_file = fopen(filename,"w");
//	fprintf(position_file,"%d\n\n",N_g);
//	for (int i = 0; i < N_g ; i++)
//	{
//		if (guests[i].type == 0)
//			fprintf(position_file,parameters.adsorbate[0].c_str());
//		if (guests[i].type == 1)
//			fprintf(position_file,parameters.adsorbate[1].c_str());
//		fprintf(position_file," %f %f %f\n",guests[i].x,guests[i].y,guests[i].z);
//	}
//}
//
////// adaptive algorithm for the translation step size
////inline void adapt_delta(double & ratio_old, double & delta_old, int & N_attempted_translations, int & N_accepted_translations, Parameters & parameters, int block_no)
////{
////	double ratio_new = (double) N_accepted_translations / (double) N_attempted_translations;
////	cout << "Block " << block_no << endl;
////	cout << "ratio new = " << ratio_new << endl;
////	if (ratio_new > 0.0)
////	{
////		if (block_no == 0) // get a starting point
////		{
////			if (ratio_new >= 0.2) // if too many translations accepted, make it bigger
////				parameters.delta = parameters.delta * 10.0;
////			else // if too few accepted, make smaller
////				parameters.delta = parameters.delta / 10.0;
////		}
////		else
////		{
////			double slope = (ratio_new - ratio_old)/(parameters.delta - delta_old); // parameters.delta is delta_new
////			parameters.delta += (0.2 - ratio_new) / slope; // update delta
////		}
////	}
////	else
////	{
////		parameters.delta = parameters.delta/ 10.0;
////	}
////	// reset count of accetpances
////	N_attempted_translations = 0; N_accepted_translations = 0;
////	// store new as old
////	ratio_old = ratio_new;
////	delta_old = parameters.delta;
////}
//
//void build_adsorbate_LJ_matrix(Forcefield forcefield, double adsorbate_sigma_matrix[2][2], double adsorbate_epsilon_matrix[2][2],Parameters parameters, ofstream & outputfile)
//{
//	assert(parameters.N_adsorbates < 3);
//	// set sigma, epsilon matrix
//	adsorbate_sigma_matrix[0][0] = parameters.sigma_guest[0]; adsorbate_sigma_matrix[1][1] = parameters.sigma_guest[1]; adsorbate_sigma_matrix[0][1] = (parameters.sigma_guest[0] + parameters.sigma_guest[1])/2; adsorbate_sigma_matrix[1][0] = (parameters.sigma_guest[0] + parameters.sigma_guest[1])/2;
//	adsorbate_epsilon_matrix[0][0] = parameters.epsilon_guest[0]; adsorbate_epsilon_matrix[1][1] = parameters.epsilon_guest[1]; adsorbate_epsilon_matrix[0][1] = sqrt(parameters.epsilon_guest[0] * parameters.epsilon_guest[1]); adsorbate_epsilon_matrix[1][0] = sqrt(parameters.epsilon_guest[0] * parameters.epsilon_guest[1]);
//	// print to outputfile
//	outputfile << "\tNumber of adsorbates = " << parameters.N_adsorbates << endl;
//	for (int n_a = 0 ; n_a < parameters.N_adsorbates; n_a ++)
//	{
//		if (n_a == 0)
//			outputfile << "\tAdsorbate " << n_a << ": " << parameters.adsorbate[0] << endl;
//		if (n_a == 1)
//			outputfile << "\tAdsorbate " << n_a << ": " << parameters.adsorbate[1] << endl;
//	}
//	outputfile << "Adosrbate Sigma matrix: ";
//	for (int i = 0 ; i < parameters.N_adsorbates; i++)
//	{
//		outputfile << "\n";
//		for (int j = 0 ; j < parameters.N_adsorbates; j++)
//		{
//			outputfile<< adsorbate_sigma_matrix[i][j] << "\t";
//		}
//	}
//	outputfile << "\nAdsorbate epsilon matrix: ";
//	for (int i=0;i<parameters.N_adsorbates;i++)
//	{
//		outputfile << "\n";
//		for (int j=0;j<parameters.N_adsorbates;j++)
//		{
//			outputfile << adsorbate_epsilon_matrix[i][j] << "\t";
//		}
//	}
//	outputfile << "\n";
//}

int main(int argc, char *argv[])
{
    if (! ((argc == 4) | (argc == 6))) {
        printf("Run as ./gcmc $structure $adsorbate0 $fugacity0(Pa) $adsorbate1 $fugactiy1(Pa)\nAdsorbate1 stuff is optional\n");
        exit(EXIT_FAILURE);
    }
    
    GCMCParameters parameters;
    parameters.frameworkname = argv[1];
    parameters.adsorbate[0] = argv[2];
    parameters.fugacity[0] = atof(argv[3]);
    if (argc == 5) {
        parameters.adsorbate[1] = argv[4];
        parameters.fugacity[1] = atof(argv[5]);
    }

    readsimulationinputfile(parameters);
    if (parameters.verbose) printf("Read simulation.input\n");
    triple_int uc_reps = readunitcellreplicationfile(parameters.frameworkname);
    parameters.replication_factor_a = uc_reps.arg1;
    parameters.replication_factor_b = uc_reps.arg2;
    parameters.replication_factor_c = uc_reps.arg3;
    if (parameters.verbose) printf("Read .uc replication file\n");
    
    parameters.adsorbateMW[0] = get_adsorbate_MW(parameters.adsorbate[0]);
    parameters.adsorbateMW[1] = get_adsorbate_MW(parameters.adsorbate[1]);

    //
    // Construct forcefield and framework objects
    //
    Forcefield forcefield(parameters.forcefieldname);
    if (parameters.verbose) printf("Constructed Forcefield object\n");
    Framework framework(parameters.frameworkname);
    parameters.N_framework_atoms = framework.noatoms;
    if (parameters.verbose) printf("Constructed Framework object\n");

//    // grab sigma/epsilon of adsorbate
//    get_guest_FF_params_from_Forcefield(forcefield, parameters.adsorbate, parameters);
//    if (parameters.verbose) printf("Fetched adsorbate FF parameters\n");
//	Parameters parameters; 
//	parameters.framework_name = argv[1];
//	parameters.fugacity[0] = atof(argv[2]); //Pa
//	if (argc == 4) parameters.fugacity[1] = atof(argv[3]); //Pa
//	// SET UP OUTPUT FILE and GRID FILE
//	ofstream outputfile;
//	char outputfilename[512]="output_files/";
//	strcat(outputfilename,parameters.framework_name.c_str());
//	outputfile.open(std::strcat(outputfilename,"_gcmc_output.txt"));
//	write_time(outputfile);
//
//	Framework framework("../sim_data/cssrs/"+parameters.framework_name+".cssr",outputfile); // construct framework object **tested.
//	Forcefield forcefield; // data structure for force field
//	import_sim_settings(parameters, outputfile, framework , parameters.framework_name, forcefield,false); // import the data, print stuff to outptufile
//	// build adsorbate LJ matrix
//	double adsorbate_sigma_matrix[2][2]; double adsorbate_epsilon_matrix[2][2];
//	build_adsorbate_LJ_matrix(forcefield, adsorbate_sigma_matrix, adsorbate_epsilon_matrix,parameters,outputfile);
//	outputfile << "\tDelta step for GCMC (A) = " << parameters.delta << endl;
//	outputfile << "\tSample every # cycles: " << parameters.sample_every << endl;
//	outputfile << "\tEquilibrium trials: " << parameters.equilibrium_trials << endl;
//	outputfile << "\tprob exchange: " << parameters.p_exchange << endl;
//	outputfile << "\tprob move: " << parameters.p_move << endl;
//	if (parameters.p_identity_change > 0.0)
//	{
//		outputfile << "\tprob ID change: " << parameters.p_identity_change << endl;
//		if (parameters.N_adsorbates == 1)
//		{
//			cerr << "Dude, why are you doing an ID change for single component?" << endl;
//			outputfile << "Dude, why are you doing an ID change for single component?" << endl;
//			exit(EXIT_FAILURE);
//		}
//	}
//  assert (parameters.p_exchange + parameters.p_move + parameters.p_identity_change < 1.0000001);
//  assert (parameters.p_exchange + parameters.p_move + parameters.p_identity_change > 0.9999999);
//
//	//
//	// IMPORT ENERGY GRIDS
//	// + pocket blocking if enabled.
//	//
//	int * N_x = (int *) calloc(2,sizeof(int)); int * N_y = (int *) calloc(2,sizeof(int));	int * N_z = (int *) calloc(2,sizeof(int));
//	double * dx_f = (double *) calloc(2,sizeof(double)); double * dy_f = (double *) calloc(2,sizeof(double));	double * dz_f = (double *) calloc(2,sizeof(double));
//	double * energy_grid0; double * energy_grid1;
//	bool pocket_block_verbose_mode = false;
//	for (int n_c = 0 ; n_c < parameters.N_adsorbates ; n_c ++)
//	{
//		if (parameters.debugmode) cout << "Importing energy grid " << n_c << endl;
//		char gridfilename[512]="../sim_data/grids/";
//		strcat(gridfilename,parameters.framework_name.c_str());
//		strcat(gridfilename,"_");	strcat(gridfilename,parameters.adsorbate[n_c].c_str());
//		strcat(gridfilename,"_"); strcat(gridfilename,parameters.forcefield.c_str());
//
//		ifstream gridfile(strcat(gridfilename,".txt")); // input file stream
//		if (gridfile.fail())
//		{
//			outputfile << "File " << gridfilename << "could not be imported." <<endl;
//			cerr << "File " << gridfilename << " could not be imported." <<endl;
//			exit(EXIT_FAILURE);
//		}
//		string line;
//		getline(gridfile,line);
//		istringstream this_line(line);
//		this_line >> N_x[n_c] >> N_y[n_c] >> N_z[n_c];
//		// create grid import format
//		if (n_c == 0)
//			energy_grid0 = (double *) calloc(N_x[0] * N_y[0] * N_z[0],sizeof(double));
//		if (n_c == 1)
//			energy_grid1 = (double *) calloc(N_x[1] * N_y[1] * N_z[1],sizeof(double));
//
//		int i=0; int j=0;
//		while(getline(gridfile,line))
//		{
//			istringstream this_line(line);
//			for (int k=0;k<N_z[n_c];k++) // each line is a pencil of z's
//			{
//				int index_here = k + j*N_z[n_c] + i*N_y[n_c]*N_z[n_c];
//				if (n_c == 0)
//					this_line >> energy_grid0[index_here];
//				if (n_c == 1)
//					this_line >> energy_grid1[index_here];
//			}
//			j+=1; // each line is a pencil of z's for a particular y... so update y index
//			if (j==N_y[n_c])
//			{
//				i += 1; // at the end of the z-y sheet, we go to a new x.
//				j = 0; // j goes back to zero
//			}
//		}
//		assert(i==N_x[n_c]); // assert that we are at the end of the file
//		assert(j==0);
//		gridfile.close();
//		// grid spacing in fractional coordinates  (how it was done in writegrid)
//		dx_f[n_c] = 1.0/(N_x[n_c]-1); // grid spacings
//		dy_f[n_c] = 1.0/(N_y[n_c]-1);
//		dz_f[n_c] = 1.0/(N_z[n_c]-1);
//		outputfile << "Size of energy grid " << n_c << " : " << N_x[n_c]*N_y[n_c]*N_z[n_c] * sizeof(double) / (1024.0 * 1024.0) << " MB\n";
//		outputfile << "\tNumber of grid points in energy grid " << N_x[n_c] << " by " << N_y[n_c] << " by " << N_z[n_c] <<endl;
//		outputfile << "\tFractional grid spacing: " << dx_f[n_c] << " by " << dy_f[n_c] << " by " << dz_f[n_c] <<endl;
//		if (parameters.debugmode) printf("energy grid %d imported successfully.\n",n_c);
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
//////// OLD WAY.
////////	char gridfilename[256]="grids/";
////////	strcat(gridfilename,parameters.framework_name.c_str());
////////	strcat(gridfilename,"_");	strcat(gridfilename,parameters.adsorbate.c_str());
////////	strcat(gridfilename,"_"); strcat(gridfilename,parameters.forcefield.c_str());
////////	ifstream gridfile(strcat(gridfilename,".txt")); // input file stream
////////	if (gridfile.fail())
////////	{
////////		outputfile << "File " << gridfilename << "could not be imported." <<endl;
////////		cerr << "File " << gridfilename << " could not be imported." <<endl;
////////		exit(EXIT_FAILURE);
////////	}
////////	string line;
////////	getline(gridfile,line);
////////	istringstream this_line(line);
////////	this_line >> N_x[0] >> N_y[0] >> N_z[0];
////////		dx_f[0] = 1.0/(N_x[0]-1); // grid spacings
////////		dy_f[0] = 1.0/(N_y[0]-1);
////////		dz_f[0] = 1.0/(N_z[0]-1);
////////	// create grid import format
////////	double * energy_grid0 = (double *) malloc(N_x[0] * N_y[0] * N_z[0] *sizeof(double));
////////	double * energy_grid1 = (double *) malloc(N_x[0] * N_y[0] * N_z[0] *sizeof(double));
////////
////////	int i=0; int j=0;
////////	while(getline(gridfile,line))
////////	{
////////		istringstream this_line(line);
////////		for (int k=0;k<N_z[0];k++) // each line is a pencil of z's
////////		{
////////			int index_here = k + j*N_z[0] + i*N_y[0]*N_z[0];
////////			this_line >> energy_grid0[index_here];
////////		}
////////		j+=1; // each line is a pencil of z's for a particular y... so update y index
////////		if (j==N_y[0])
////////		{
////////			i += 1; // at the end of the z-y sheet, we go to a new x.
////////			j = 0; // j goes back to zero
////////		}
////////	}
////////	assert(i==N_x[0]); // assert that we are at the end of the file
////////	assert(j==0);
////////	outputfile << "energy grid imported successfully." <<endl;
////
////		//
////		//  POCKET BLOCKING: block inaccessible pockets by overwriting values in * energy_grid with high numbers in blocked regions
////		//
////		if (parameters.pocket_blocking)
////		{
////			outputfile << "Pocket blocking enabled with temperature " << parameters.T << " K." << endl;
////			/* --- RICH'S CODE BEGINS HERE --- */
////			double time_before = read_timer();
////			if (pocket_block_verbose_mode == true)
////			{
////				cout << "Pocket blocking beginning. Write a cube for before" << endl;
////				outputfile << "VERBOSE MODE FOR POCKETBLOCKING!" <<endl;
////				if (n_c == 0)
////					write_cube("before_blocking_0", framework, &parameters, energy_grid0, N_x[n_c],N_y[n_c],N_z[n_c]); //just so you can see the grid before ...
////				if (n_c == 1)
////					write_cube("before_blocking_1", framework, &parameters, energy_grid1, N_x[n_c],N_y[n_c],N_z[n_c]); //just so you can see the grid before ...
////				cout << "Cube written." << endl;
////			}
////			int num_threads = 8; //set this wherever you like or make it a parameter
////			int num_pockets = -1;
////			if (n_c == 0)
////				num_pockets = find_and_block_pockets(framework, energy_grid0, parameters.T, num_threads, false, N_x[n_c],N_y[n_c],N_z[n_c], parameters);
////			if (n_c == 1)
////				num_pockets = find_and_block_pockets(framework, energy_grid1, parameters.T, num_threads, false, N_x[n_c],N_y[n_c],N_z[n_c], parameters);
////			outputfile << parameters.adsorbate[n_c] << " : found " << num_pockets << " inaccessible pockets" << endl;
////			// energy_grid is passed as a poitner, so it is modified now if accessible pockets were there.
////			if (pocket_block_verbose_mode == true)
////			{
////				cout << "Pocket blocking finished. Write a cube for after" << endl;
////				double time_after = read_timer();
////				outputfile << "Time spent to find and block pockets using " << num_threads << " threads: " << time_after-time_before << endl;
////				if (n_c == 0)
////					write_cube("after_blocking", framework, &parameters, energy_grid0, N_x[n_c],N_y[n_c],N_z[n_c]); // ... and after pocket blocking
////				if (n_c == 1)
////					write_cube("after_blocking", framework, &parameters, energy_grid1, N_x[n_c],N_y[n_c],N_z[n_c]); // ... and after pocket blocking
////				cout << "Cube written." << endl;
////			}
////			/* --- RICH'S CODE ENDS HERE --- */
////		}
////	}
//
//  // store framework positions in fractional coordinates, sigma and epsilon via mixing rule. replicate framework if needed
//	particle_f * framework_atoms = (particle_f *) malloc(parameters.N_framework_atoms * sizeof(particle_f));
//	for (int f = 0; f < framework.noatoms; f++)
//	{
//		// keep framework_atoms on device as fractional coordinates [0,1]^3
//		framework_atoms[f].x = framework.atoms[f].x;
//		framework_atoms[f].y = framework.atoms[f].y;
//		framework_atoms[f].z = framework.atoms[f].z;
//
//		string ff_label = framework.atoms[f].label;
//		bool found_frameworkatom=false;
//		double ep; double sig; // for storing framework epsilon and sigma
//		for (int q = 0; q < forcefield.no_interactions; q++)
//		{
//			if (forcefield.label[q] == ff_label)
//			{
//				ep = forcefield.epsilon[q];
//				sig = forcefield.sigma[q];
//				found_frameworkatom = true;
//				break;
//			}
//		}
//		if (found_frameworkatom == false)
//		{
//			outputfile << "Error: Could not find framework atom " << ff_label << " in the forcefield file." <<endl;
//			exit(0);
//		}
////		if (forcefield.mixing_rule=="Lorentz-Berthelot")
////		{
//		// mixing rules
//		framework_atoms[f].sigma = (parameters.sigma_guest[0] + sig)/2.0;
//		framework_atoms[f].epsilon = sqrt(parameters.epsilon_guest[0] * ep);
////		}
////		else
////		{
////			outputfile << "Error: Don't have " << forcefield.mixing_rule << " mixing rule."<<endl;
////			exit(0);
////		}
//	}
////	print_interaction_parameters(framework,framework_atoms,outputfile,parameters);
//
//	if (parameters.debugmode)
//		cout << "Imported framework atoms\n";
//
//	//
//	// CHECK ENERGY GRID WITH QUICK CALCULATION IN SERIAL CODE
//	//
//	double energy_serial_code = 0.0; // energy at test point calcuated from serial code
//	// test at point:
//	// 0.458379 0.782972 0.0868405 1.044343,yf=0.324289,zf=0.945021
//	double x_f_check = .044343; double y_f_check = 0.324289; double z_f_check = 0.945021;
//	getEnergyAtPoint(x_f_check,y_f_check,z_f_check,framework_atoms,parameters,energy_serial_code);
//	double energy_grid_interpolation = E_gf(x_f_check, y_f_check, z_f_check,N_x[0], N_y[0], N_z[0],dx_f[0], dy_f[0], dz_f[0],energy_grid0);
//	outputfile << "CHECKING GRID at test point x_f,y_f,z_f "<< x_f_check<< " , " << y_f_check << " , " << z_f_check <<endl;
//	outputfile << "Energy calculated from serial code: " << energy_serial_code <<endl;
//	outputfile << "Energy calculated from grid interpolation: " << energy_grid_interpolation <<endl;
//	x_f_check = 0.5; y_f_check = 0.8; z_f_check = 0.4;
//	getEnergyAtPoint(x_f_check,y_f_check,z_f_check,framework_atoms,parameters,energy_serial_code);
//	energy_grid_interpolation = E_gf(x_f_check, y_f_check, z_f_check,N_x[0], N_y[0], N_z[0],dx_f[0], dy_f[0], dz_f[0],energy_grid0);
//	outputfile << "CHECKING GRID at test point x_f,y_f,z_f "<< x_f_check<< " , " << y_f_check << " , " << z_f_check <<endl;
//	outputfile << "Energy calculated from serial code: " << energy_serial_code <<endl;
//	outputfile << "Energy calculated from grid interpolation: " << energy_grid_interpolation <<endl;
//
//	//
//	// GCMC stats and guest positions (frac and cartesian)
//	//
//	GCMC_stats statistics;
//	double volume = parameters.unitcell_volume * parameters.replication_factor_a * parameters.replication_factor_b * parameters.replication_factor_c;
//	// guests array has both components (Xe and Kr).  type indicates identity
//	particle_g * guests = (particle_g *) malloc(MAX_GUESTS * sizeof(particle_g));
//  outputfile << "\tSize of guests = " << MAX_GUESTS * sizeof(particle_g) / (1024.0 * 1024.0) << " MB\n";
//
//	//
//	// SET UP RANDOM NUMBER GENERATOR
//	//
//	  // uniform real in [0,1]
//	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
////	std::default_random_engine generator (seed);
//	std::mt19937 generator (seed);
//	std::uniform_real_distribution<double> uniform01 (0.0,1.0); // random nos in 0,1
//	std::uniform_int_distribution<int> uniformint(0,10); // initialize with bogus range, will change later.
//	std::uniform_int_distribution<int> type_generator(0,parameters.N_adsorbates - 1); // initialize with bogus range, will change later.
//
////	std::cout << "some random numbers between 0.0 and 1.0: " << std::endl;
////	for (int i=0; i<10; ++i)
////		std::cout << uniform01(generator) << std::endl;
//		// uniform int in [0,z]
//
//	//
//	//  RUN GCMC
//	//
//	double start_of_sim_time=read_timer();
//	outputfile << "Starting simulation" << endl;
//	// initialize stuff
//	int * N_g = (int *) calloc(2,sizeof(int)); // current number of guests
//	N_g[0] = 0; N_g[1] = 0;
//	statistics.N_samples = 0;
//	statistics.N_g_avg[0] = 0.0;	statistics.N_g_avg[1] = 0.0;
//	statistics.N_moves = 0;	statistics.N_insertions = 0; 	statistics.N_deletions = 0;
//	statistics.N_move_trials = 0;	statistics.N_insertion_trials = 0; 	statistics.N_deletion_trials = 0;
//	statistics.guest_guest_energy_avg = 0.0;
//	statistics.framework_guest_energy_avg= 0.0;
//
//
//	// initialize
//	double E_gg_this_cycle = 0.0; double E_gf_this_cycle = 0.0; // assumes we start with zero particles
//	int cycle_counter = 0;
//	int adsorbate_index_list[2][MAX_GUESTS] = { -1 }; // keep indices of particles here
//	int N_g_total = 0;
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
//				statistics.N_insertion_trials += 1;
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
//				if (parameters.debugmode)
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
//					if (parameters.debugmode) cout << "\tInsertion accepted. " << endl;
//					statistics.N_insertions += 1;
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
//				if (parameters.debugmode) cout << "Deletion Trial." << endl;
//				statistics.N_deletion_trials += 1;
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
//						statistics.N_deletions += 1;
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
//						if (parameters.debugmode)
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
//				if (parameters.debugmode) cout << "Translation Trial." << endl;
//				statistics.N_move_trials += 1;
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
//						statistics.N_moves += 1; N_accepted_translations += 1;
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
//					statistics.N_ID_swap_trials += 1;
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
//						statistics.N_ID_swaps += 1;
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
//				statistics.N_samples += 1;
//				statistics.N_g_avg[0] += N_g[0]; statistics.N_g_avg[1] += N_g[1];
//				statistics.N_g2_avg[0] += N_g[0] * N_g[0]; statistics.N_g2_avg[1] += N_g[1] * N_g[1];
//				statistics.guest_guest_energy_avg += E_gg_this_cycle;
//				statistics.framework_guest_energy_avg += E_gf_this_cycle;
//			}
//			if (parameters.debugmode == true)
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
//	statistics.guest_guest_energy_avg = 1.0 * statistics.guest_guest_energy_avg / statistics.N_samples;
//	statistics.framework_guest_energy_avg = 1.0 * statistics.framework_guest_energy_avg / statistics.N_samples;
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
//	outputfile << "N_insertions: " << statistics.N_insertions << " / " << statistics.N_insertion_trials << endl;
//	outputfile << "N_deletions: " << statistics.N_deletions << " / " << statistics.N_deletion_trials << endl;
//	outputfile << "N_moves: " << statistics.N_moves << " / " << statistics.N_move_trials << endl;
//	outputfile << "N_ID_Swaps: " << statistics.N_ID_swaps << " / " << statistics.N_ID_swap_trials << endl;
//	outputfile << "N_samples: " << statistics.N_samples << endl;
//
//	outputfile << "\n\nRESULTS\n=======\n";
//	outputfile << "\tTemperature (K): " << parameters.T << endl;
//	for (int n_c = 0 ; n_c < parameters.N_adsorbates ; n_c ++)
//	{
//		string adsorbate;
//		if (n_c == 0)
//			adsorbate = parameters.adsorbate[0];
//		if (n_c == 1)
//			adsorbate = parameters.adsorbate[1];
//		statistics.N_g_avg[n_c] = 1.0 * statistics.N_g_avg[n_c] / statistics.N_samples;
//		statistics.N_g2_avg[n_c] = 1.0 * statistics.N_g2_avg[n_c] / statistics.N_samples;
//		double N_confidence_bound = sqrt(statistics.N_g2_avg[n_c] - statistics.N_g_avg[n_c] * statistics.N_g_avg[n_c])/sqrt(1.0 * statistics.N_samples); // sigma / sqrt(N)
//		outputfile << "Adsorbate " << n_c << ": " << adsorbate <<endl;
//		outputfile << "\tFugacity (Pa) = " << parameters.fugacity[n_c] << endl;
//		outputfile << "\t<N_g> (" << adsorbate << ") = " << 1.0 * statistics.N_g_avg[n_c] / parameters.replication_factor_a / parameters.replication_factor_b / parameters.replication_factor_c << " +/- " << N_confidence_bound << " molecules per unit cell" << endl;
//		outputfile << "\t<N_g> (" << adsorbate << ", moles/m3) = " << statistics.N_g_avg[n_c]/volume/6.022e-7 << endl;
//		outputfile << "\t<N_g> (" << adsorbate << ", moles/kg framework) = " << statistics.N_g_avg[n_c]/volume/6.022e-7/framework.density << endl;
//	}
//	outputfile << "<E_gg> (K) = " << statistics.guest_guest_energy_avg <<endl; //* 8.314 / 1000.0<< endl;
//	outputfile << "\t(kJ/mol) = " << statistics.guest_guest_energy_avg * 8.314 / 1000.0<< endl;
//	outputfile << "<E_gf> (K) = " << statistics.framework_guest_energy_avg <<endl; //* 8.314 / 1000.0<< endl;
//	outputfile << "\t(kJ/mol) = " << statistics.framework_guest_energy_avg * 8.314 / 1000.0<< endl;
}
