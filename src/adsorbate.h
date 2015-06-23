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
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class Adsorbate {
    public:
        int type;  // adsorbate identity
        int nbeads;  // number of beads in adsorbate
        matrix<double> bead_xyz; // Cartesian coordinates of beads stored column-wise (at most 2 at this point)
        int beadtypes[2];  // consists of max 2 beads (leave second blank if one bead...)

        Adsorbate()
            :bead_xyz(3, 2)  // preallocate size as 3 by 2 (at most 2 beads)
        {
            // Constructor (bogus defaults)
            type = -1;
            nbeads = 0;
        }

        void translate_by_Cartesian_vector(double dx, double dy, double dz) {
            // translate by Cartesian coords (dx, dy, dz)
            for (int b = 0; b < nbeads; b++) {
                bead_xyz(0, b) += dx;
                bead_xyz(1, b) += dy;
                bead_xyz(2, b) += dz;
            }
        }

        void print_info() {
            printf("\nAdsorbate type %d.\n", type);
            printf("\tNumber of beads = %d\n", nbeads);
            for (int b = 0; b < nbeads; b++)
               printf("\tBead %d. Type %d. (%f, %f, %f)\n", b, beadtypes[b], bead_xyz(0, b), bead_xyz(1, b), bead_xyz(2, b));
            printf("\n");
        }
};

void PerformUniformRandomRotation(Adsorbate & adsorbate) {
    // rotate adsorbate molecule
    double R[3][3];  // rotation matrix
}

std::map<std::string, int> GetUniqueBeads(std::vector<std::string> adsorbatelist, bool verbose) {
    // return map of unique beads
    std::map<std::string, int> uniquebeads;
    int nuniquebeads = 0;
    
    for (int a = 0; a < adsorbatelist.size(); a++) {
        char adsorbatefilename[1024];
        sprintf(adsorbatefilename, "../data/adsorbates/%s.adsorbate", adsorbatelist[a].c_str());
        std::ifstream adsorbatefile(adsorbatefilename);
        if (adsorbatefile.fail()) {
            printf("adsorbate info file not found: data/adsorbates/%s.adsorbate\n", adsorbatelist[a].c_str());
            exit(EXIT_FAILURE);
        }
        // Get number of beads in this adsorbate
        int nbeads;
        std::string line;
        std::getline(adsorbatefile, line);
        std::istringstream linestream(line);
        linestream >> nbeads;
        if (verbose)
            printf("Adsorbate %s: %d beads.\n", adsorbatelist[a].c_str(), nbeads);

        // Get the bead names
        std::getline(adsorbatefile, line);
        linestream.str(line); linestream.clear();
        for (int b = 0; b < nbeads; b++) {
            std::string beadname;
            linestream >> beadname;
            if (verbose)
                printf("\t%s\n", beadname.c_str());
            if (uniquebeads.find(beadname) == uniquebeads.end()) {
                // if bead not found, add it
                uniquebeads[beadname] = nuniquebeads;
                nuniquebeads += 1;
            }
        }

        adsorbatefile.close();
    }

    if (verbose) {
        printf("%d unique beads.\n", nuniquebeads);

        for (std::map<std::string, int>::iterator it=uniquebeads.begin(); it!=uniquebeads.end(); ++it)
            std::cout << it->first << " => " << it->second << '\n';
    }
    return uniquebeads;
};

std::vector<Adsorbate> GetAdsorbateTemplates(std::vector<std::string> adsorbatelist, std::map<std::string, int> uniquebeads, bool verbose) {
    // Get adsorbate templates for adsorbate molecules by reading files
    // Ensure first bead is at 0, 0, 0
    // return vector of adsorbate templates
    std::vector<Adsorbate> adsorbatetemplates(adsorbatelist.size());

    for (int a = 0; a < adsorbatelist.size(); a++) {
        adsorbatetemplates[a].type = a;

        char adsorbatefilename[1024];
        sprintf(adsorbatefilename, "../data/adsorbates/%s.adsorbate", adsorbatelist[a].c_str());
        std::ifstream adsorbatefile(adsorbatefilename);
        if (adsorbatefile.fail()) {
            printf("adsorbate info file not found: data/adsorbates/%s.adsorbate\n", adsorbatelist[a].c_str());
            exit(EXIT_FAILURE);
        }
        // Get number of beads in this adsorbate
        int nbeads;
        std::string line;
        std::getline(adsorbatefile, line);
        std::istringstream linestream(line);
        linestream >> nbeads;
        adsorbatetemplates[a].nbeads = nbeads; 
        if (verbose)
            printf("Adsorbate %s: %d beads.\n", adsorbatelist[a].c_str(), nbeads);

        std::vector<std::string> beadnames;
       
        // Get the bead names
        std::getline(adsorbatefile, line);
        linestream.str(line); linestream.clear();
        for (int b = 0; b < nbeads; b++) {
            std::string beadname;
            linestream >> beadname;
            if (verbose)
                printf("\t%s\n", beadname.c_str());
            beadnames.push_back(beadname);
        }
         
        // get postitions
        getline(adsorbatefile, line);  // "bead positions.csv"
        getline(adsorbatefile, line);  // "xyz"
        for (int b = 0; b < nbeads; b++) {
            std::getline(adsorbatefile, line);
            std::cout << line << std::endl;
            linestream.str(line); linestream.clear();
            
            adsorbatetemplates[a].beadtypes[b] = uniquebeads[beadnames[b]];

            linestream >> adsorbatetemplates[a].bead_xyz(0, b); 
            linestream >> adsorbatetemplates[a].bead_xyz(1, b);
            linestream >> adsorbatetemplates[a].bead_xyz(2, b); 
        }

        adsorbatefile.close();
    }
    return adsorbatetemplates;
}

#endif /* ADSORBATE_H_ */
