/*
Class file for methods used in more classes that are directly related to the rrt rather then being general whatever usable everywhere like in the support_methods class.
*/
#include "world.h"
#include "obstacle.h"

extern const int CHECK_ONLY_BLOCKING_SPHERES; 
extern const int CHECK_WITH_BLOCKING_SPHERES; 
extern const int DONT_CHECK_WITH_BLOCKING_SPHERES;

bool is_in_obstacle(double* loc_coord, double radius, int check_with_blocking_spheres);
bool is_in_obstacle(std::array<double, 3> loc_coord, double radius, int check_with_blocking_spheres);
bool is_in_obstacle_custom_frame(double* loc_coord, double radius, int check_with_blocking_spheres, int frame);
//bool is_in_obstacle(array<double, 3> loc_coord, double radius, int check_with_blocking_spheres);
double approximate_radius(double* coordinates, double default_radius, double precision, double default_approximation_step);
double approximate_radius_custom_frame(double* coordinates, double default_radius, double precision, double default_approximation_step, int frame);
void cut_subtree_in_main_tree(shared_ptr<vertex> first_delet_node, bool has_debugging_output, bool is_nearest_neighbor_search_structure_rebuilt);
void rebuild_kd_tree(bool is_initialized, bool copy_stuff, bool only_valid_in_curframe);
int add_to_tree(double* loc_coord, int tree_index, MPNN::MultiANN<double>* kdTree);
void build_tree(int is_global);
bool test_path_noncolliding_static(shared_ptr<Path> tested_path);

int print_colliding(std::map<int ,shared_ptr<vertex>>& map);
void purge_nodes_in_blocking_spheres(shared_ptr<Path> optimizedTunnel);

void set_path_bottleneck(shared_ptr<Path> path);
void set_path_length(shared_ptr<Path> path);