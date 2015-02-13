#include<assert.h>
#include<string>
#include<cstdlib> // for "exit" and "malloc"
#include "datatypes.h"
#include "Framework.h"
#include<cstring>
#include <stdio.h>
#include<vector>

inline int GetEnergyGridIndex(int i, int j, int k, GridInfo grid_info) {
    // returns index in energy_grid corresponding to gridpoint index i,j,k
    return k + j * grid_info.N_z + i * grid_info.N_y * grid_info.N_z;
}

struct Grid_point_info {
    bool status_found;
    int segment_id; // = -1  if inaccessible, 0,1,2,3,... for other segments.
};

struct Connection {// for segment connectivity
    int from;
    int to;
    int direction;
};

struct Unit_cell_id { // to keep track of which unit cell we are in.
    int i, j, k;
};

struct Grid_id { // store a grid point
    int i, j, k; 
    int global_id;
};

// define directions
int X_DIR = 1; int Y_DIR = 2; int Z_DIR = 3; // positive for forward

// makes channel_found true if loop was found; false if connections exhausted without finding loop.
void recursively_navigate_connections(int segment_no, 
                        std::vector<Unit_cell_id> & unit_cell_each_segment, 
                        bool * segment_visit_status, 
                        bool * edge_visit_status, 
                        std::vector<Connection> & connections, 
                        int N_segments, 
                        int N_connections, 
                        bool * channel_found)
{
    // we are exploring segement number segment_no.
    segment_visit_status[segment_no] = true; // now this has been visited.
    // loop over all UNVISTED connections
    for (int c = 0; c < N_connections; c ++) {
        if (edge_visit_status[c] == false) {
            // self-connection with this segment, kk this is a channel!
            if ((connections[c].to == segment_no) && (connections[c].from == segment_no)) {
                edge_visit_status[c] = true;
                (* channel_found) = true;
            }
            // connection between segment_no and a different segment.
            else if ((connections[c].to == segment_no) || (connections[c].from == segment_no)) {
                edge_visit_status[c] = true;
                // we want to travel from THIS segment to the other, so reverse if needed.
                int other_segment = connections[c].to;
                int direction = connections[c].direction;
                if (connections[c].to == segment_no) {
                    other_segment = connections[c].from;
                    direction = -1 * direction; // reverse direction
                }
                // kk so going from segment segment_no to other_segment in the direction direction.
                // update unit cell we will be in
                Unit_cell_id this_uc = unit_cell_each_segment.at(segment_no);
                Unit_cell_id new_uc = this_uc;
                if (direction == X_DIR)
                    new_uc.i  ++;
                if (direction == Y_DIR)
                    new_uc.j  ++;
                if (direction == Z_DIR)
                    new_uc.k  ++;
                if (direction == -X_DIR)
                    new_uc.i  --;
                if (direction == -Y_DIR)
                    new_uc.j  --;
                if (direction == -Z_DIR)
                    new_uc.k  --;
                // if segment not visited, walk along a different segment
                if ( segment_visit_status[other_segment] == false) {
                    unit_cell_each_segment.at(other_segment) = new_uc;
                    recursively_navigate_connections(other_segment, unit_cell_each_segment, segment_visit_status, edge_visit_status, connections, N_segments, N_connections, channel_found);
                }
                // if it is has been visited already, don't visit it again, but do check whether the unit cell is the same as the old one
                else {
                    // get unit cell of other segment
                    Unit_cell_id old_uc = unit_cell_each_segment.at(other_segment);
                    // if we arrive at the same segment in a different unit cell, this is a channel!
                    if ((new_uc.i != old_uc.i) || (new_uc.j != old_uc.j) || (new_uc.k != old_uc.k))
                        (*channel_found) = true;
                }
            }
        } // end if edge has not been visited
    } // end loop over all connections
}

// we observed a connection. should we add it to the connections list, or does it already exist?
void note_new_connection(int from_segment, int to_segment, int direction, std::vector<Connection> & connections) {
    bool should_push = true; // don't push a redundant connection to the connections std::vector...
    for (int i = 0; i < connections.size(); i++) {
    // look at all present connections to see if it exists already
        if ((to_segment == connections.at(i).to) && (from_segment == connections.at(i).from)) {
            // does a connection between this node already exist?
            if ((to_segment == from_segment) || (direction == connections.at(i).direction)) { 
                // if to and from are same, direction does not matter. if direction from x to y same, we hv this connection...
                should_push = false;
                break;
            }
        }
    }
    if (should_push) {
        // if we should push this connection
        Connection c;
        c.to = to_segment;
        c.from = from_segment;
        c.direction = direction;
        connections.push_back(c); // add connection
    }
}

int FindAndBlockPockets(double * energy_grid, GridInfo grid_info, double temperature, GCMCParameters parameters)
{
    if (parameters.verbose) printf("Starting blocking pockets routine.\n");
    int grid_size = grid_info.numtotalpts; // number of grid points
    double energy_threshold = 15.0 * temperature; // energy < energy_threshold ==> inaccessible
    if (parameters.verbose) printf("Energy threshold = %f\n",energy_threshold);

    Grid_point_info * grid_pt_info = (Grid_point_info *) malloc(grid_size * sizeof(Grid_point_info)); // status grid. stores status of each grid point. + segment to which it belongs
    //
    // initialize all points as FAR, unless inaccessible. Then: segment_id = 0.
    //
    for (int i = 0; i < grid_size; i++) {
        if (energy_grid[i] >= energy_threshold) {
            // grid point is inaccessible
            grid_pt_info[i].status_found = true;
            grid_pt_info[i].segment_id = -1; // segment -1 is for inaccessible
        }
        else
            grid_pt_info[i].status_found = false;
    }
    if (parameters.verbose) printf("Allocated grid point status\nStarting flood fill..\n");

    int N_segments = 0; // count segments, excludes inacessible segment (0)
    //
    // FLOOD FILL loop over all grid points to assign segments. ignore periodic bcs for now
    //
    for (int i = 0; i < grid_info.N_x; i++) {
        for (int j = 0; j < grid_info.N_y; j++) {
            for (int k = 0; k < grid_info.N_z; k++) {
                // declare grid point
                Grid_id grid_pt; grid_pt.i = i; grid_pt.j = j; grid_pt.k = k; grid_pt.global_id = GetEnergyGridIndex(i,j,k,grid_info);
                if (grid_pt_info[grid_pt.global_id].status_found == false) {
                    N_segments += 1; // we have arrived at a new segment!
                    if (parameters.debugmode) printf("\tFound a new segment\nN_segments = %d\n",N_segments);
                    std::vector<Grid_id> stack_ids; // stack with IDs unknown
                    stack_ids.push_back(grid_pt);
                    while (! stack_ids.empty()) {
                        // get another grid point from the stack
                        Grid_id another_grid_pt = stack_ids.at(stack_ids.size() - 1); // get an element from the stack.
                        stack_ids.pop_back(); // remove last element
                        // this grid pt is adjacent to the grid_pt, thus belongs to same segment.
                        grid_pt_info[another_grid_pt.global_id].status_found = true;
                        grid_pt_info[another_grid_pt.global_id].segment_id = N_segments - 1; // assign segment number
                        // push each adjacent point to the stack, as this adjacent point must belong to the same segment!
                        // x, left and right
                        if (another_grid_pt.i != 0) {
                            // declare pt
                            Grid_id adjacent_pt;
                            adjacent_pt.i = another_grid_pt.i - 1;
                            adjacent_pt.j = another_grid_pt.j;
                            adjacent_pt.k = another_grid_pt.k;
                            adjacent_pt.global_id = GetEnergyGridIndex(adjacent_pt.i,adjacent_pt.j,adjacent_pt.k,grid_info);
                            if (grid_pt_info[adjacent_pt.global_id].status_found == false) {
                                // push it to the stack
                                stack_ids.push_back(adjacent_pt);
                            }
                        }
                        if (another_grid_pt.i != grid_info.N_x - 1) {
                            // declare pt
                            Grid_id adjacent_pt;
                            adjacent_pt.i = another_grid_pt.i + 1;
                            adjacent_pt.j = another_grid_pt.j;
                            adjacent_pt.k = another_grid_pt.k;
                            adjacent_pt.global_id = GetEnergyGridIndex(adjacent_pt.i,adjacent_pt.j,adjacent_pt.k,grid_info);
                            if (grid_pt_info[adjacent_pt.global_id].status_found == false) {
                                // push it to the stack
                                stack_ids.push_back(adjacent_pt);
                            }
                        }
                        // y , left and right
                        if (another_grid_pt.j != 0) {
                            // declare pt
                            Grid_id adjacent_pt;
                            adjacent_pt.i = another_grid_pt.i;
                            adjacent_pt.j = another_grid_pt.j - 1;
                            adjacent_pt.k = another_grid_pt.k;
                            adjacent_pt.global_id = GetEnergyGridIndex(adjacent_pt.i,adjacent_pt.j,adjacent_pt.k,grid_info);
                            if (grid_pt_info[adjacent_pt.global_id].status_found == false)
                            {
                                // push it to the stack
                                stack_ids.push_back(adjacent_pt);
                            }
                        }
                        if (another_grid_pt.j != grid_info.N_y - 1) {
                            // declare pt
                            Grid_id adjacent_pt;
                            adjacent_pt.i = another_grid_pt.i;
                            adjacent_pt.j = another_grid_pt.j + 1;
                            adjacent_pt.k = another_grid_pt.k;
                            adjacent_pt.global_id = GetEnergyGridIndex(adjacent_pt.i,adjacent_pt.j,adjacent_pt.k,grid_info);
                            if (grid_pt_info[adjacent_pt.global_id].status_found == false)
                            {
                                // push it to the stack
                                stack_ids.push_back(adjacent_pt);
                            }
                        }
                        // z , left and right
                        if (another_grid_pt.k != 0) {
                            // declare pt
                            Grid_id adjacent_pt;
                            adjacent_pt.i = another_grid_pt.i;
                            adjacent_pt.j = another_grid_pt.j;
                            adjacent_pt.k = another_grid_pt.k - 1;
                            adjacent_pt.global_id = GetEnergyGridIndex(adjacent_pt.i,adjacent_pt.j,adjacent_pt.k,grid_info);
                            if (grid_pt_info[adjacent_pt.global_id].status_found == false)
                            {
                                // push it to the stack
                                stack_ids.push_back(adjacent_pt);
                            }
                        }
                        if (another_grid_pt.k != grid_info.N_z - 1) {
                            // declare pt
                            Grid_id adjacent_pt;
                            adjacent_pt.i = another_grid_pt.i;
                            adjacent_pt.j = another_grid_pt.j;
                            adjacent_pt.k = another_grid_pt.k + 1;
                            adjacent_pt.global_id = GetEnergyGridIndex(adjacent_pt.i,adjacent_pt.j,adjacent_pt.k,grid_info);
                            if (grid_pt_info[adjacent_pt.global_id].status_found == false)
                            {
                                // push it to the stack
                                stack_ids.push_back(adjacent_pt);
                            }
                        }
                    } // end while stack size != 0
                } // end IF GRID POINT IS FAR
    }}} // end loop over grid points

    // at this point, all points should be identified. Now, need to implement periodic BC's and check adjacent ones.
    if (parameters.verbose) printf("Flood fill done\n");
    if (parameters.verbose) printf("There are %d segments present.\n", N_segments);
    if (parameters.verbose) printf("Does it matter if grid includes redundant point i.e. 0 and 1 fractioanl coords?\n");

    if (parameters.verbose) printf("Checking that all grid points have been found a status\n");
    for (int i = 0 ; i < grid_size; i++) {
        if (grid_pt_info[i].status_found == false) {
            printf("Dammit, missed status of grid point %d\n",i);
            exit(EXIT_FAILURE);
        }
    }
    if (parameters.verbose) printf("\tYAY all grid points hv a status now\n");

    //
    // BUILD GRAPH THAT SHOWS HOW SEGMENTS ARE CONNECTED VIA PBC
    //
    // store graph and the direction of their connections
    if (parameters.verbose) printf("Extracting graph representation of segments\n");
    //
    std::vector<Connection> connections; // store connections across periodic boundary here
    //
    // loop over all faces of unit cell
    // store only positive connections (I think...)
    int i = 0; // bottom and top faces
    for (int j = 0; j < grid_info.N_y; j++) {
        for (int k = 0; k < grid_info.N_z; k++) {
            // get this point on the face and its segment
            int this_global_id = GetEnergyGridIndex(i,j,k,grid_info);
            int this_segment = grid_pt_info[this_global_id].segment_id;
            // get neighbor and its segment
            int neighbor_global_id = GetEnergyGridIndex(grid_info.N_x-1,j,k,grid_info);
            int neighbor_segment = grid_pt_info[neighbor_global_id].segment_id;
            if ((this_segment != -1) && (neighbor_segment != -1)) {
//                  note_new_connection(this_segment, neighbor_segment, x_dir_neg, connections);
                note_new_connection(neighbor_segment, this_segment, X_DIR, connections); // goes in positive direction
            }
        }
    }
    int j = 0; // front and back faces
    for (int i = 0; i < grid_info.N_x; i++) {
        for (int k = 0; k < grid_info.N_z; k++) {
            // get this point on the face and its segment
            int this_global_id = GetEnergyGridIndex(i,j,k,grid_info);
            int this_segment = grid_pt_info[this_global_id].segment_id;
            // get neighbor and its segment
            int neighbor_global_id = GetEnergyGridIndex(i,grid_info.N_y-1,k,grid_info);
            int neighbor_segment = grid_pt_info[neighbor_global_id].segment_id;
            if ((this_segment != -1) && (neighbor_segment != -1)) {
//                  note_new_connection(this_segment, neighbor_segment, y_dir_neg, connections);
                note_new_connection(neighbor_segment, this_segment, Y_DIR, connections);
            }
        }
    }
    int k = 0; // front and back faces
    for (int i = 0; i < grid_info.N_x; i++) {
        for (int j = 0; j < grid_info.N_y; j++) {
            // get this point on the face and its segment
            int this_global_id = GetEnergyGridIndex(i,j,k,grid_info);
            int this_segment = grid_pt_info[this_global_id].segment_id;
            // get neighbor and its segment
            int neighbor_global_id = GetEnergyGridIndex(i,j,grid_info.N_z-1,grid_info);
            int neighbor_segment = grid_pt_info[neighbor_global_id].segment_id;
            if ((this_segment != -1) && (neighbor_segment != -1)) {
//                  note_new_connection(this_segment, neighbor_segment, z_dir_neg, connections);
                note_new_connection(neighbor_segment, this_segment, Z_DIR, connections);
            }
        }
    } // graph has been built, characterized in terms of connections std::vector.
    int N_connections = connections.size();
    if (parameters.verbose) printf("Connections graph created. Number of connections = %d.\n",N_connections);
    for (int i = 0; i < N_connections; i++)
        if (parameters.debugmode) printf("Connection from seg %d to %d, direction = %d.\n",connections.at(i).from,connections.at(i).to,connections.at(i).direction);
    //
    // LABEL SEGMENT AS POCKET OR CHANNEL BASED ON GRAPH APPROACH
    //
    // store status of each segment
    int UNKNOWN = 1; int POCKET = 2; int CHANNEL = 3;
    if (parameters.debugmode) printf("Status IDs: UNKNOWN = %d, POCKET = %d, CHANNEL = %d\n",UNKNOWN,POCKET,CHANNEL);
    std::vector<int> segment_identification;
    std::vector<int> segments_that_are_pockets;
    for (int i = 0; i < N_segments; i++)
        segment_identification.push_back(UNKNOWN);
    int N_pockets = 0;
    // create base of unit cells for each segment
    Unit_cell_id unit_cell_BASE; unit_cell_BASE.i = 0; unit_cell_BASE.j = 0; unit_cell_BASE.k = 0;
    std::vector<Unit_cell_id> unit_cells_each_segment_BASE;
    for (int orz = 0; orz < N_segments; orz++)
        unit_cells_each_segment_BASE.push_back(unit_cell_BASE);
    // keep track of edges and segments that were visited
    bool * segment_visit_status = (bool *) malloc(N_segments * sizeof(bool));
    bool * edge_visit_status = (bool *) malloc(N_connections * sizeof(bool));
    for (int seg = 0; seg < N_segments ; seg++) {
        // navigate through connections until all are exhausted or a loop is detected.
        // set visit statues as false
        for (int count = 0; count < N_segments; count ++)
            segment_visit_status[count] = false;
        for (int count = 0; count < N_connections; count ++)
            edge_visit_status[count] = false;
        bool channel_found = false; // pass by reference to make it a pointer
        std::vector<Unit_cell_id> unit_cells_each_segment = unit_cells_each_segment_BASE; // make them all (0,0,0)
        recursively_navigate_connections(seg, unit_cells_each_segment, segment_visit_status, edge_visit_status, connections, N_segments, N_connections, & channel_found);
        if (channel_found) {
            if (parameters.verbose) printf("Found that segment %d is a channel.\n",seg);
            segment_identification.at(seg) = CHANNEL;
        }
        if (! channel_found) {
            segment_identification.at(seg) = POCKET;
            if (parameters.verbose) printf("Didn't find a channel for segment %d.\n",seg);
            segments_that_are_pockets.push_back(seg);
            N_pockets++;
        }
    }
    assert(N_pockets == segments_that_are_pockets.size());
    if (parameters.verbose) printf("Number of pockets identified = %d\n",N_pockets);

    //
    // FOR EACH POCKET, OVERWRITE ENERGY WITH HIGH VALUE
    //
    if (parameters.verbose) printf("Overwriting pockets with large energy value\n");
    double high_energy = 1e13;
    for (int pp = 0; pp < N_pockets; pp++) {
        if (parameters.verbose) printf("\tOverwriting pocket %d\n",segments_that_are_pockets.at(pp));
        // loop through all grid points. if this belongs to the pockets, overwrite
        for (int gg = 0; gg < grid_size; gg++)
        {
            if (grid_pt_info[gg].segment_id == segments_that_are_pockets.at(pp))
                energy_grid[gg] = high_energy;
        }
    }
    //
    // FREE MEMORY (IMPORTANT HERE)
    //
    if (parameters.debugmode) printf("Freeing memory\n");
    free(segment_visit_status);
    free(edge_visit_status);
    free(grid_pt_info);
    return N_pockets;
}


void WriteCube(std::string cube_name, Framework framework, GCMCParameters parameters, double * energy_grid, GridInfo grid_info) {
    // write energy grid pointer malloc'ed array to cube file
    FILE * gridfile;
    char gridfilename[512];
    sprintf(gridfilename, "data/grids/%s.cube", cube_name.c_str());
    gridfile = fopen(gridfilename, "w");
    fprintf(gridfile, "\nThis is a grid file.\n");
    fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n", 0, 0.0, 0.0, 0.0); //enforce zero atoms
    fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n", grid_info.N_x, parameters.t_matrix[0][0]/(grid_info.N_x - 1),0.0,0.0);
    fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n", grid_info.N_y, parameters.t_matrix[0][1]/(grid_info.N_y - 1), parameters.t_matrix[1][1] / (grid_info.N_y - 1), 0.0);
    fprintf(gridfile, "%d % 13.6lf % 13.6lf % 13.6lf\n", grid_info.N_z, parameters.t_matrix[0][2]/(grid_info.N_z - 1), parameters.t_matrix[1][2] / (grid_info.N_z - 1), parameters.t_matrix[2][2] / (grid_info.N_z - 1));
    for(int i = 0; i < grid_info.N_x; i++) {
        for(int j = 0; j < grid_info.N_y; j++) { // loop over y's
            int count = 0;
            for(int k = 0; k < grid_info.N_z; k++) {
                fprintf(gridfile, "% 13.6E ", energy_grid[GetEnergyGridIndex(i, j, k, grid_info)]);
                count ++;
                if( count == 6) {
                    fprintf(gridfile,"\n");
                    count=0;
                }
            }
            fprintf(gridfile, "\n"); //new line after z over
        }
    }
    fclose(gridfile);
}
