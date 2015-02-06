/*
 * pocket_blocking.h
 *   Find and block pockets from energy grid
 *  Created on: Feb 4, 2015
 *      Author: corymsimon
 */
#ifndef POCKETBLOCKING_H
#define POCKETBLOCKING_H

int find_and_block_pockets(double * energy_grid, Grid_info grid_info, double temperature, GCMCParameters parameters);

void write_cube(string cube_name, Framework framework, GCMCParameters parameters, double * energy_grid, Grid_info grid_info);
#endif
