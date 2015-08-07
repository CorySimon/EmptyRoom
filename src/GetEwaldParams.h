// for loading Ewald parameters structure
#ifndef Ewald_H
#define Ewald_H
#include<string>
#include<vector>
#include "Framework.h"
#include "datatypes.h"

std::vector<double> cross_product(std::vector<double> a, std::vector<double> b) {
    // return cross product of vectors a and b (3D only)
    std::vector<double> result(3);
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
}

double inner_product(std::vector<double> a, std::vector<double> b) {
    // dot product
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

EWaldParameters GetEwaldParams(Framework framework, bool verbose=false) {
    // Fill Ewald Parameters
    //
    EWaldParameters params; 

    //  vacuum permittivity  eps0 = 8.854187817e-12 C^2/(J-m)
    //  1 m = 1e10 A, 1 e = 1.602176e-19 C, kb = 1.3806488e-23 J/K
    //  eps0 = 8.854187817e-12 C^2/(J-m) [1 m / 1e10 A] [1 e / 1.602176e-19 C]^2 [kb = 1.3806488e-23 J/K]
    params.eps0 = 4.7622424954949676e-7;  // \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)
    
    // in Gaussian used for splitting short and long range interactions (see Berend's book)
    params.alpha = 0.025;  

    // short-range cutoff
    params.cutoff_squared = 21.0 * 21.0;

    // replications in fourier space for long-range interactions
    params.kx = 7;
    params.ky = 7;
    params.kz = 6;

    //
    // Get Reciprocal lattice vectors
    //
    // Unit cell lattice vectors from transformation matrix
    std::vector< std::vector<double> > a;
    for (int i = 0; i < 3; i++) {
        std::vector<double> ai(3);
        ai[0] = framework.t_matrix[0][i];
        ai[1] = framework.t_matrix[1][i];
        ai[2] = framework.t_matrix[2][i];
        a.push_back(ai);
    }
    
    // cross products
    std::vector<double> a1_a2 = cross_product(a[1], a[2]);
    std::vector<double> a2_a0 = cross_product(a[2], a[0]);
    std::vector<double> a0_a1 = cross_product(a[0], a[1]);
    
    // Reciprocal lattice space vectors
    std::vector< std::vector<double> > b;
    for (int i = 0; i < 3; i++) {
        std::vector<double> bi(3);
        // multiply by 2\pi here so we don't hv to do it later
        bi[0] = 2.0 * M_PI * a1_a2[i] / inner_product(a[0], a1_a2);
        bi[1] = 2.0 * M_PI * a2_a0[i] / inner_product(a[1], a2_a0);
        bi[2] = 2.0 * M_PI * a0_a1[i] / inner_product(a[2], a0_a1);
        b.push_back(bi);
        params.b1[i] = bi[0];
        params.b2[i] = bi[1];
        params.b3[i] = bi[2];
    }
    
    if (verbose) {
        printf("Primitive lattice vectors:\n");
        for (int i = 0; i < 3; i++) 
            printf("a_%d = (%f, %f, %f)\n", i, a[i][0], a[i][1], a[i][2]);
        printf("Reciprocal lattice vectors:\n");
        for (int i = 0; i < 3; i++) 
            printf("b_%d = (%f, %f, %f)\n", i, b[i][0], b[i][1], b[i][2]);
        printf("Checking orthogonality (2 pi = %f).\n\tA * B = \n", 2.0 * M_PI);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                printf("dot(a_%d, b_%d) = %f\n", i, j,
                        a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]);
            }
        }
    }
    return params;
}

#endif /* EWald_H_ */
