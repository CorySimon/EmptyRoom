#include<iostream>
#include "Framework.h"
#include<fstream>
#include<string>
#include<random>
#include<chrono>
#include<cmath>
#include<map>
#include<set>
#include <limits>
#include <complex>
#include <cstring>
#include<math.h>
#include<assert.h>
#include<cstdlib> // for "exit"
#include<sstream> // string stream
#include<vector>
#include "adsorbate.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
namespace boo = boost::numeric::ublas;

struct bs {
    int a;
    std::map<int, std::string> gah;
};

int main()
{   
   Adsorbate xe(1);
//   xe.print_info();
//   xe.nbeads = 1;
//   xe.type = 0;
//   xe.bead_xyz(0,0) = 1.0;
//   xe.bead_xyz(1,0) = 2.0;
//   xe.bead_xyz(2,0) = 3.0;
//   xe.beadtypes[0] = 0;
//   xe.print_info();
//   printf("matrix %zu by %zu\n", xe.bead_xyz.size1(), xe.bead_xyz.size2());
//
//
//   printf("translate by (1,2,3)\n");
//   xe.translate_by_Cartesian_vector(1,2,3);
//   xe.print_info();
//   std::vector<std::multiset<int>> adsorbateIDs_by_type;
//   adsorbateIDs_by_type.resize(2);
//   adsorbateIDs_by_type[0].insert(10);
//   adsorbateIDs_by_type[0].insert(20);
//   adsorbateIDs_by_type[0].insert(15);
//   std::multiset<int>::iterator it = adsorbateIDs_by_type[0].begin();
//   std::cout << *it <<std::endl;
//
//   
//   boo::matrix<double> m(3,3);
//   m(0,1)=5;
//   boo::matrix_column<boo::matrix<double> > mc(m, 1);
//
//   boo::vector<double> x(3);
//   x[0] =.2; x[1] = .4; x[2] =3.;
//   boo::vector<double> y(3);
//   y[0] =.2; y[1] = .5; y[2] =3.;
//   boo::vector<double> z= x+ mc;
//   for (int i=0;i<3;i++)
//       printf("z(%d)=%f\n", i, z[i]);
//
   std::vector<std::string> adsorbatelist(1);
   adsorbatelist[0] = "Xe";
//   adsorbatelist[1] = "CH2CH2";

   std::map<std::string, int> uniquebeads = GetBeadMap(adsorbatelist, true);

   std::vector<Adsorbate> adsorbatetemplates = GetAdsorbateTemplates(adsorbatelist, uniquebeads, true);
   for (int i = 0; i < adsorbatelist.size(); i++) {
       adsorbatetemplates[i].print_info();
   }
//  
////  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
////    std::mt19937 generator(seed);
////    std::normal_distribution<double> distribution(0.0,1.0);
////    boo::vector<double> r = GetUniformVectorOnSphere(generator, distribution);
////    printf("r = (%f, %f, %f)\n", r[0], r[1], r[2]);
////
////    for (int rot = 0; rot< 10; rot++) {
////        printf("rotation %d\n", rot);    
////        boo::vector<double> r = GetUniformVectorOnSphere(generator, distribution);
////        printf("norm of vector on sphere = %f\n", boo::norm_2(r));
////        PerformUniformRandomRotation(adsorbatetemplates[1], generator, distribution);
////        adsorbatetemplates[1].print_info();
////        double dx = adsorbatetemplates[1].bead_xyz(0,0) - adsorbatetemplates[1].bead_xyz(0,1);
////        double dy = adsorbatetemplates[1].bead_xyz(1,0) - adsorbatetemplates[1].bead_xyz(1,1);
////        double dz = adsorbatetemplates[1].bead_xyz(2,0) - adsorbatetemplates[1].bead_xyz(2,1);
////        printf("Bond length = %f\n", sqrt(dx*dx+dy*dy+dz*dz));
////    }
    return 0;
}
