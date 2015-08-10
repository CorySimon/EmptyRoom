/*
 * adsorbate.h
 *
 *  Created on: Feb 2, 2015
 *      Author: corymsimon
 */

#ifndef ADSORBATE_H_
#define ADSORBATE_H_
#include<string>
#include<vector>
#include<random>
#include <boost/numeric/ublas/matrix.hpp>
// #include <boost/geometry/geometry.hpp>
namespace boo = boost::numeric::ublas;


bool OutsideUnitCell(double & x_f, double & y_f, double & z_f,
                     std::vector<int> & uc_reps) {
    // check if particle is outside bounding box
    if ((x_f > uc_reps[0]) | 
        (y_f > uc_reps[1]) | 
        (z_f > uc_reps[2]) | 
        (x_f < 0.0) | (y_f < 0.0) | (z_f < 0.0)) {
        return true;
    }
    else
        return false;
}

class Adsorbate {
    public:
        int type;  // adsorbate identity

        //
        // Lennard Jones beads
        //
        int nbeads;  // number of LJ beads in adsorbate
        boo::matrix<double> bead_xyz; // Cartesian coordinates of beads stored column-wise
        boo::matrix<double> bead_xyz_f; // Fractional coordinates of beads stored column-wise
        std::vector<int> beadtypes;  // Integer ID of bead type

        //
        // Point charges
        //
        int ncharges;
        std::vector<double> charges;  // Charge on bead
        boo::matrix<double> charge_xyz; // Cartesian coordinates of charges stored column-wise
        boo::matrix<double> charge_xyz_f; // Fractional coordinates of charges stored column-wise

        //
        // Constructor. Sets sizes of matrices/vectors
        //
        Adsorbate(int nbeads_, int ncharges_)
            :bead_xyz(3, nbeads_),  // preallocate size as 3 by nbeads
             bead_xyz_f(3, nbeads_),  // preallocate size as 3 by nbeads
             charge_xyz(3, ncharges_),
             charge_xyz_f(3, ncharges_),
             charges(ncharges_),
             beadtypes(nbeads_)
        {
            nbeads = nbeads_;
            ncharges = ncharges_;
            // give bogus types for safety
            type = -1;
            for (int i = 0; i < nbeads_; i++)
                beadtypes[i] = -1;
        }

        void translate_by_Cartesian_vector(boo::vector<double> & dx, 
                                           boo::matrix<double> & t_matrix,
                                           boo::matrix<double> & inv_t_matrix,
                                           std::vector<int> & uc_reps) {
            // translate by Cartesian vector dx
            
            //
            // update Cartesian coords of LJ beads and charges
            //
            for (int b = 0; b < nbeads; b++) {
                bead_xyz(0, b) += dx[0];
                bead_xyz(1, b) += dx[1];
                bead_xyz(2, b) += dx[2];
            }
            for (int c = 0; c < ncharges; c++) {
                charge_xyz(0, c) += dx[0];
                charge_xyz(1, c) += dx[1];
                charge_xyz(2, c) += dx[2];
            }
            
            //
            // update fractional coords
            //
            bead_xyz_f = boo::prod(inv_t_matrix, bead_xyz);
            charge_xyz_f = boo::prod(inv_t_matrix, charge_xyz);

            //
            //  Apply periodic BCs
            //
            // if FIRST LJ bead outside simulation box, move ENTIRE guest
            bool outside_simbox = false;
            // Loop over fractional coords. If outside [0, uc_reps[i]]...
            for (int i = 0; i < 3; i++) {
                if (bead_xyz_f(i, 0) > 1.0 * uc_reps[i]) {
                    for (int b = 0; b < nbeads; b++)
                        bead_xyz_f(i, b) -= 1.0 * uc_reps[i];
                    for (int c = 0; c < ncharges; c++)
                        charge_xyz_f(i, c) -= 1.0 * uc_reps[i];
                    outside_simbox = true;
                }
                else if (bead_xyz_f(i, 0) < 0.0) {
                    for (int b = 0; b < nbeads; b++)
                        bead_xyz_f(i, b) += 1.0 * uc_reps[i];
                    for (int c = 0; c < ncharges; c++)
                        charge_xyz_f(i, c) += 1.0 * uc_reps[i];
                    outside_simbox = true;
                }
            }
            if (outside_simbox) {
                // update Cartesian
                bead_xyz = boo::prod(t_matrix, bead_xyz_f);
                charge_xyz = boo::prod(t_matrix, charge_xyz_f);
            }
        }

        void print_info() {
            printf("\nAdsorbate type %d.\n", type);
            printf("\t%d Lennard-Jones beads\n", nbeads);
            for (int b = 0; b < nbeads; b++)
               printf("\tBead %d. Type %d. x=(%f, %f, %f)\n", b, beadtypes[b], bead_xyz(0, b), bead_xyz(1, b), bead_xyz(2, b));
            printf("\t%d Point charges\n", ncharges);
            for (int c = 0; c < ncharges; c++)
                printf("\tCharge %d. Charge = %f, x= (%f, %f, %f)\n", c, charges[c], charge_xyz(0, c), charge_xyz(1, c), charge_xyz(2, c));
            printf("\n");
        }
};

boo::vector<double> GetUniformVectorOnSphere(std::mt19937 & generator, 
                                        std::normal_distribution<double> & std_normal_distn) {
    // returns uniform random vector on sphere
    // Generates 3 normally distributed numbers with zero mean
    // Normalizes result
    boo::vector<double> r(3);
    double norm = 0.0;  // L2 norm
    while (norm < 0.00001) {
        // don't let norm get too small
        r(0) = std_normal_distn(generator);
        r(1) = std_normal_distn(generator);
        r(2) = std_normal_distn(generator);

        norm = boo::norm_2(r);
    }

    r = r / norm;
    
    return r;
}

boo::vector<double> cross_product(boo::vector<double> & a, boo::vector<double> & b) {
    // return cross product of vectors a and b (3D only)
    boo::vector<double> result(3);
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
}

void PerformUniformRandomRotation(Adsorbate & adsorbate, 
                                  std::mt19937 & generator, 
                                  std::normal_distribution<double> & std_normal_distn) {
    // rotate adsorbate molecule. WARNING: only updates Cartesian coords, not fractional
    // MUST have first bead at (0, 0, 0) to preserve bond lengths
    if ((adsorbate.bead_xyz(0, 0) != 0.0) | (adsorbate.bead_xyz(1, 0) != 0.0) | (adsorbate.bead_xyz(1, 0) != 0.0)) {
        // TODO: remove later
        printf("First bead must be at origin.\n");
        adsorbate.print_info();
        exit(EXIT_FAILURE);
    }
    
    //
    // Get rotation matrix
    //
    boo::matrix<double> R(3, 3);  // rotation matrix
    // first column is unform vector on sphere
    boo::vector<double> r1 = GetUniformVectorOnSphere(generator, std_normal_distn);
    R(0, 0) = r1(0); 
    R(1, 0) = r1(1); 
    R(2, 0) = r1(2);
    // second column is orthogonal
    double norm = 0.0;
    boo::vector<double> r2;
    while (norm < 0.00001) {
        r2 = GetUniformVectorOnSphere(generator, std_normal_distn);
        r2 = r2 - boo::inner_prod(r2, r1) * r1;  // subtract off component along r1
        norm = boo::norm_2(r2);
    }
    r2 = r2 / norm;
    R(0, 1) = r2(0); 
    R(1, 1) = r2(1); 
    R(2, 1) = r2(2);
    // third column is cross product of two to give another orthogonal vector
    boo::vector<double> r3 = cross_product(r1, r2);
    r3 = r3 / boo::norm_2(r3);
    R(0, 2) = r3(0); 
    R(1, 2) = r3(1); 
    R(2, 2) = r3(2);
//    printf("r1 dot r2 = %f\n", boo::inner_prod(r2, r1));
//    printf("r1 dot r3 = %f\n", boo::inner_prod(r1, r3));
//    printf("r2 dot r3 = %f\n", boo::inner_prod(r2, r3));
//    for (int i = 0; i<3;i++)
//        printf("%f %f %f\n", R(i,0), R(i,1), R(i,2));
    
    //
    // Rotate Cartesian coords
    //
    adsorbate.bead_xyz = boo::prod(R, adsorbate.bead_xyz);
    adsorbate.charge_xyz = boo::prod(R, adsorbate.charge_xyz);
}

std::map<std::string, int> GetBeadMap(std::vector<std::string> adsorbatelist, bool verbose) {
    // return map of unique beads. bead label --> integer mapping
    std::map<std::string, int> beadlabel_to_int;
    int nuniquebeads = 0;
    
    for (int a = 0; a < adsorbatelist.size(); a++) {
        char adsorbatefilename[1024];
        sprintf(adsorbatefilename, "data/adsorbates/%s.adsorbate", adsorbatelist[a].c_str());
        std::ifstream adsorbatefile(adsorbatefilename);
        if (adsorbatefile.fail()) {
            printf("adsorbate info file not found...: %s\n", adsorbatefilename);
            exit(EXIT_FAILURE);
        }

        //
        // Get number of beads in this adsorbate
        //
        int nbeads;
        std::string line;
        std::getline(adsorbatefile, line);
        std::istringstream linestream(line);
        linestream >> nbeads;
        if (verbose)
            printf("Adsorbate %s: %d beads.\n", adsorbatelist[a].c_str(), nbeads);
        if (!(nbeads > 0)) {
            printf("Number of beads not greater than zero. Something wrong in format of adsorbate info file.\n");
            exit(EXIT_FAILURE);
        }

        // Get the bead names
        std::getline(adsorbatefile, line);  // x point charges
        std::getline(adsorbatefile, line);  // ---------------
        std::getline(adsorbatefile, line);  // Lennard-Jones beads:
        std::getline(adsorbatefile, line);  // Name x y z
        for (int b = 0; b < nbeads; b++) {
            std::getline(adsorbatefile, line);  // [Name x y z]_{bead i}
            linestream.str(line); linestream.clear();
            std::string beadname;
            linestream >> beadname;
            if (verbose)
                printf("\t%s\n", beadname.c_str());
            if (beadlabel_to_int.find(beadname) == beadlabel_to_int.end()) {
                // if bead not found, add it
                beadlabel_to_int[beadname] = nuniquebeads;
                nuniquebeads += 1;
            }
        }

        adsorbatefile.close();
    }

    if (verbose) {
        printf("%d unique Lennard-Jones beads.\n", nuniquebeads);

        for (std::map<std::string, int>::iterator it=beadlabel_to_int.begin(); it!=beadlabel_to_int.end(); ++it)
            std::cout << it->first << " => " << it->second << '\n';
    }
    return beadlabel_to_int;
};

std::vector<Adsorbate> GetAdsorbateTemplates(std::vector<std::string> adsorbatelist, std::map<std::string, int> beadlabel_to_int, bool verbose) {
    // Get adsorbate templates for adsorbate molecules by reading files
    // Ensure first bead is at 0, 0, 0
    // return vector of adsorbate templates
    int numadsorbates = adsorbatelist.size();
    std::vector<Adsorbate> adsorbatetemplates;

    for (int a = 0; a < adsorbatelist.size(); a++) {

        char adsorbatefilename[1024];
        sprintf(adsorbatefilename, "data/adsorbates/%s.adsorbate", adsorbatelist[a].c_str());
        std::ifstream adsorbatefile(adsorbatefilename);
        if (adsorbatefile.fail()) {
            printf("adsorbate info file not found: %s\n", adsorbatefilename);
            exit(EXIT_FAILURE);
        }
        
        //
        // Get number of LJ beads and point charges in this adsorbate
        //
        int nbeads, ncharges;
        // LJ beads
        std::string line;
        std::getline(adsorbatefile, line);
        std::istringstream linestream(line);
        linestream >> nbeads;
        // point charges 
        std::getline(adsorbatefile, line);
        linestream.str(line); linestream.clear();
        linestream >> ncharges;
        if (verbose)
            printf("Adsorbate %s: %d LJ beads, %d charges.\n", adsorbatelist[a].c_str(), nbeads, ncharges);
        
        //
        // construct adsorbate
        //
        Adsorbate adsorbate(nbeads, ncharges);
        adsorbate.type = a;
        
        //
        // Get the bead names and positions
        //
        std::string beadname;
        std::getline(adsorbatefile, line);  // -----
        std::getline(adsorbatefile, line);  // Lennard-Jones beads:
        std::getline(adsorbatefile, line);  // Name x y z
        for (int b = 0; b < nbeads; b++) {
            std::getline(adsorbatefile, line);
            linestream.str(line); linestream.clear();

            linestream >> beadname >> adsorbate.bead_xyz(0, b) >> adsorbate.bead_xyz(1, b) >> adsorbate.bead_xyz(2, b);
            adsorbate.beadtypes[b] = beadlabel_to_int[beadname];
        }

        //
        // Get point charges
        //
        std::getline(adsorbatefile, line);  // -----
        std::getline(adsorbatefile, line);  // Point charges:
        std::getline(adsorbatefile, line);  // x y z charge
        for (int c = 0; c < ncharges; c++) {
            std::getline(adsorbatefile, line);
            linestream.str(line); linestream.clear();
            linestream >> adsorbate.charge_xyz(0, c) >> adsorbate.charge_xyz(1, c) >> adsorbate.charge_xyz(2, c) >> adsorbate.charges[c];
        }
        
        //
        // Ensure first LJ bead is at origin
        //
        for (int b = 0; b < nbeads; b++) {
            adsorbate.bead_xyz(0, b) -= adsorbate.bead_xyz(0, 0);
            adsorbate.bead_xyz(1, b) -= adsorbate.bead_xyz(1, 0);
            adsorbate.bead_xyz(2, b) -= adsorbate.bead_xyz(2, 0);
        }
        for (int c = 0; c < ncharges; c++) {
            adsorbate.charge_xyz(0, c) -= adsorbate.bead_xyz(0, 0);
            adsorbate.charge_xyz(1, c) -= adsorbate.bead_xyz(1, 0);
            adsorbate.charge_xyz(2, c) -= adsorbate.bead_xyz(2, 0);
        }
        
        // 
        // Check for charge neutrality
        //
        double total_charge = 0.0;
        for (int c = 0; c < ncharges; c++)
            total_charge += adsorbate.charges[c];
        
        if ((total_charge > 0.0001) | (total_charge < -0.0001)) {
            printf("Net charge of adsorbate %d is %f != 0.0\n", adsorbate.type, total_charge);
            for (int c = 0; c < ncharges; c++)
                printf("\tPoint charge %d charge: %f\n", c, adsorbate.charges[c]);
        }

        adsorbatefile.close();

        adsorbatetemplates.push_back(adsorbate);
    }

    return adsorbatetemplates;
}

#endif /* ADSORBATE_H_ */
