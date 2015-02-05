/*
 * pocket_blocking.h
 *   Find and block pockets from energy grid
 *  Created on: Feb 4, 2015
 *      Author: corymsimon
 */
#ifndef CORYFLOODFILL_H
#define CORYFLOODFILL_H
int find_and_block_pockets(double * energy_grid, int N_x, int N_y, int N_z, double temperature, GCMCParameters parameters);
void write_cube(string cube_name, Framework framework, GCMCParameters *p, double *energy_grid, int N_x, int N_y, int N_z);
#endif