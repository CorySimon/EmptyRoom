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

// An adsorbate consists of a set of beads
struct Bead {
    int type;  // identity of bead
    
    // Cartesian coords
    double x;
    double y;
    double z;
};

class Adsorbate {
    public:
        int type;  // adsorbate identity
        int nbeads;  // number of beads in adsorbate
        Bead beads[2];  // consists of max 2 beads (leave second blank if one bead...)
        double centerofmass[3];  // center of mass

        Adsorbate() {
            // Constructor (bogus defaults)
            type = -1;
            nbeads = -1;
        }

        void translate_by_Cartesian_vector(double dx, double dy, double dz) {
            // translate by Cartesian coords (dx, dy, dz)
            for (int b = 0; b < nbeads; b++) {
                beads[b].x += dx;
                beads[b].y += dy;
                beads[b].z += dz;
            }
        }

        void print_info() {
            printf("\nAdsorbate type %d.\n", type);
            printf("\tNumber of beads = %d\n", nbeads);
            for (int b = 0; b < nbeads; b++)
               printf("\tBead %d. Type %d. (%f, %f, %f)\n", b, beads[b].type, beads[b].x, beads[b].y, beads[b].z);
            printf("\n");
        }
};

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
        getline(adsorbatefile, line);
        std::istringstream linestream(line);
        linestream >> nbeads;
        if (verbose)
            printf("Adsorbate %s: %d beads.\n", adsorbatelist[a].c_str(), nbeads);

        // Get the bead names
        getline(adsorbatefile, line);
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

#endif /* ADSORBATE_H_ */
