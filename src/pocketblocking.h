/*
 * pocket_blocking.h
 *   Find and block pockets from energy grid
 *  Created on: Feb 4, 2015
 *      Author: corymsimon
 */
#ifndef POCKETBLOCKING_H
#define POCKETBLOCKING_H

int FindAndBlockPockets(double * energy_grid, GridInfo grid_info, double temperature, GCMCParameters parameters);

void WriteCube(std::string cube_name, Framework framework, GCMCParameters parameters, double * energy_grid, GridInfo grid_info);
#endif
