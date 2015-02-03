/*
 * Writes Potential energy grid of adsorbate inside unit cell of nanoporous material
 */
#include <stdio.h>
#include <stdlib.h>
#include<cuda.h>
using namespace std;
#include "readsettings.h"
#include "datatypes.h"
#include "Framework.h"
#include "Forcefield.h"
#include "write_to_outputfile.h"

void host_frac_to_cart(double T_matrix[][3],
		double x_f, double y_f,double z_f,
		double & x, double & y, double & z) {
	// compute Cartesian coordinates from fractional
    x = T_matrix[0][0] * x_f + T_matrix[0][1] * y_f + T_matrix[0][2] * z_f;
    y = T_matrix[1][0] * x_f + T_matrix[1][1] * y_f + T_matrix[1][2] * z_f;
    z = T_matrix[2][0] * x_f + T_matrix[2][1] * y_f + T_matrix[2][2] * z_f;
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
    parameters.adsorbateMW = get_adsorbate_MW_of(parameters.adsorbate);

    readsimulationinputfile(parameters);
    readunitcellreplicationfile(parameters, parameters.frameworkname);

    //
    // Construct forcefield and framework objects
    //
    Forcefield forcefield(parameters.forcefieldname);
    Framework framework(parameters.frameworkname);

    get_guest_FF_params_from_Forcefield(forcefield, parameters.adsorbate, parameters); // get sig/eps of adsorbate

    //
    //  Write settings to outputfile
    //
    FILE * outputfile;
    char outputfilename[512] = "output_files/";
    strcat(outputfilename, parameters.frameworkname.c_str());
    strcat(outputfilename, "_"); strcat(outputfilename, parameters.adsorbate.c_str()); strcat(outputfilename, "_grid.out");
    outputfile = fopen(outputfilename, "w");
    write_settings_to_outputfile(outputfile, parameters, framework);

}
