#ifndef TUNNEL_OPTIMIZATION
#define TUNNEL_OPTIMIZATION

#include "world.h"
#include "obstacle.h"
#include "support_methods.h"

double* center_node(shared_ptr<vertex> centered_node, 
	shared_ptr<vertex> previous_node, shared_ptr<vertex> next_node, double* coordinate_vector, double increment, bool are_neighbors_centered, bool rewrite_node, bool return_node_info);
int path_optimization(shared_ptr<Path> path);


#endif