#include <iostream>
#include <iomanip>
#include <fstream>
#include <dirent.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <vector>
#include <chrono>
#include <map>
#include <sstream>
#include <memory>

#include "world.h"
#include "obstacle.h"
#include "MPNN/include/DNN/multiann.h"
#include "MPNN/include/DNN/ANN.h"
#include "support_methods.h"
#include "path_duplication_check.h"
#include "tunnel_optimization.h"
#include "rrt_methods.h"
#include "file_management_methods.h"

using namespace std::chrono;
using namespace ozcollide;

//to do: nonholonomicky hledani, pridani kd stromu
//parameters:
//functions, structs definitions
int find_nearest_vertex(double* loc_coord, MPNN::MultiANN<double>* kdTree);
void write_paths_to_pdb(const char* filename, vector<shared_ptr<vertex>>& path);
double* collision_function(double x1, double y1, double x2, double y2, double min_step_col);
int add_to_tree(double* loc_coord, int tree_index, MPNN::MultiANN<double>* kdTree);
void build_tree(int is_global);
void print_vertex_coordinates(vertex* vert);
shared_ptr<vertex> step(shared_ptr<vertex> new_vertex, shared_ptr<vertex> parent, double goal_distance, 
  double collision_check_step, int parent_index);
void backtrack_process_path(shared_ptr<vertex> new_vertex );
void step_interpolated(shared_ptr<vertex> new_vertex, shared_ptr<vertex> nearest_neighbor_vertex, double max_step, 
  double interpolation_step, int parent_index);

double MAX_STEP_DISTANCE = 1;


bool found = false;
bool is_finished = false;


std::vector <int> duplicate_count; //keeps track of path count in "clusters"

//--- used to keep track for tunnel validity in frames
vector <int> paths_count;
vector <bool> path_found_in_frame;
//-----

double maximum_increase = 1.1; //aximum increase in blocking sphere size due to many duplicates, in percents (1.5 = 150 percent bigger sphere)
bool step_success_flag;        //used to determine succesful add to local tree, if fails global is tryied or if global doesn't exist, iteration fail

bool LOCALIZED_MODE = false; 
//kept for debugging purposes---
bool are_duplicated_pdbs_created = false;
std::map <int, shared_ptr<vertex>> blocking_spheres_in_vertices;
std::vector <std::vector<shared_ptr<Path>>> duplicated_paths;
std::vector <shared_ptr<Path>> cluster_paths;
//vector <vertex*> visualise_centring;
//---


//-----variables for finding program statistics
int nearest_vertex_search_count = 0;
int num_of_duplicates = 0;
//-----


//------program parameters  //consider using txt file
const double DUPLICATION_CHECK_CONSTANT = 1.2; //ratio of inter-tunnel distance to tunnel length, if lower, tunnel is considered a duplicate of existing tunnel
const double BLOCKING_SPHERE_RADIUS_COEFFICIENT = 1.2;
//******these are used for debugging purposes, set them to some absurd number for normal run
const int MAX_DUPLICATES = 10000;               //num of duplicates before run terminates //
const int MAX_TUNNELS = 10000;                   //maximum number of tunnels before run terminates
//******
//------

//
void add_cluster_path(shared_ptr<Path> added_path, int path_index){
  if(duplicated_paths.size() <= path_index){
    std::vector<shared_ptr<Path>> new_cluster;
    new_cluster.push_back(added_path);
    //cout << "cluster " << new_cluster.size()<<endl;
    duplicated_paths.push_back(new_cluster);
    //cout << "dupl" << duplicated_paths.size()<<endl;
    return;
  } 

  if(!path_found_in_frame[path_index]){
   duplicated_paths[path_index].push_back(added_path);
  }
}

//obvious
shared_ptr<vertex> load_parent_vertex_from_tree(shared_ptr<vertex> cur_vertex){
  return local_priority_kdTree_coordinates->get_node_pointer(cur_vertex->get_parent_index());
}

//frame is valid if both tested node and it's child are valid in it and and the tested node's is same or smaller "is previous"
//when locating parent node, we can't "travel into the future", the node should be in the same or older frame than the actual. If it is valid in more than one,obviously the most recent is used.
/*
int find_valid_frame(shared_ptr<vertex> cur, shared_ptr<vertex> child){
  int nearest_frame = child->get_frame_index();
  //cout << "nearest frame " << nearest_frame << endl;
  while(nearest_frame != -1){

    if(cur->is_valid_in_frame(nearest_frame)){
      return nearest_frame;
    }
    nearest_frame = child->find_previous_frame(nearest_frame);
    //cout << "nearest frame search " << nearest_frame << endl;
  }
  return -1;
}
*/

//sets node into already found frame, deletes all children other than child in path, returns pointer to node added to path
//deletes other children pointers so the path has only the path child, 
shared_ptr<vertex> copy_node_tree_to_path(shared_ptr<Path> path, shared_ptr<vertex> parent_tree_node, shared_ptr<vertex> child_node, int child_first_frame){


  shared_ptr<vertex> parent_cur_node = parent_tree_node->copy_without_structure_pointers();
  parent_cur_node->set_last_valid_frame(child_first_frame);
  parent_cur_node->add_child_pointer(child_node, false);
  child_node->set_parent_pointer(parent_cur_node);
  path->add_node(parent_cur_node);

  return parent_cur_node;
}

//creates a copy of vertices involved in path from found endpoint to initial point, the path must be copied in case that it's part is deleted
shared_ptr<Path> backtrack(shared_ptr<vertex> endpoint){ 
  shared_ptr<Path> path = make_shared<Path>(-1, endpoint->get_index(), get_current_frame());
  shared_ptr<vertex> cur = load_parent_vertex_from_tree(endpoint);

  //cout << "cur index " << cur->get_index() << endl;
  //cout << "cur frame " << cur->get_last_valid_frame() << endl;

  if(cur->get_last_valid_frame() != endpoint->get_last_valid_frame()) {cerr << "error in backtrack, endpoint and it's parents last frames doesn't match! " << endl; exit(0);} //these are more for testing, shoudn't happen
  if(endpoint->get_last_valid_frame() != get_current_frame()) {cerr << "error in backtrack, endpoint isn't in current frame! " << endl; exit(0);}
  if(endpoint->get_children_count() > 0) {cerr << "error in backtrack, endpoint has children! " << endl; exit(0);}

  shared_ptr<vertex> endpoint_copied = endpoint->copy_without_structure_pointers(); //
  path->add_node(endpoint_copied);

  shared_ptr<vertex> child_pointer = copy_node_tree_to_path(path, cur, endpoint_copied, endpoint_copied->get_first_frame());
  

  while(cur->get_parent_index() != -1){ //cur is previous node from the tree
    cur = load_parent_vertex_from_tree(cur);
    int child_first_frame = child_pointer->get_first_frame();
    //sanity check
    if(!cur->is_valid_in_frame(child_first_frame)){
      cout << "rrt_main_file.cpp, backtrack, insane " << endl;
      create_segfault();
    }
    child_pointer = copy_node_tree_to_path(path, cur, child_pointer, child_pointer->get_first_frame());
/*
    cout << "index " << child_pointer->get_index() << endl;
    cout << "first frame " <<  child_pointer->get_first_frame() << endl;
    cout << "last frame " <<  child_pointer->get_last_valid_frame() << endl;
    cout << "inactive " <<  child_pointer->is_inactive() << endl;
*/

  }
  path->set_beginning_index(cur->get_index()); 
  return path;
  
}

//initializes k-d tree, puts initial vertex with parent index -1, symbolizing it hasn't got a parent, into tree map
void init_start(){
  build_tree(LOCAL);
  local_priority_kdTree_coordinates = new Tree();
  double location_coordinates[3] = {qix, qiy, qiz};
  shared_ptr<vertex> start = make_shared<vertex>(location_coordinates, 0, probe_radius, get_current_frame(), get_current_frame(), true);
  local_priority_kdTree_coordinates->add_node(start);
  
  tree_index = 0;  //use as a static variable in tree::add_node;
  add_to_tree(start->get_location_coordinates(), tree_index, local_priority_kdTree);
  tree_index++; 
  found = false;
  is_local_tree_initialized = true;
}


  

//function for generating next sample
//todo - use percentage 
double* sample(int iterations_counter){
  double vertX; double vertY; double vertZ;
  double* loc_coord = new double[3];
  int counter = 0;
    while(true){
      vertX = (((double) rand() / (double) RAND_MAX) * (world_size_x)) + world_size_x_shift;
      vertY = (((double) rand() / (double) RAND_MAX) * (world_size_y)) + world_size_y_shift;
      vertZ = (((double) rand() / (double) RAND_MAX) * (world_size_z)) + world_size_z_shift;
      loc_coord[0] = vertX; loc_coord[1] =  vertY; loc_coord[2] = vertZ;
      double r = ((double) rand() / (RAND_MAX));
      counter++;
      //cout << "test " << counter << " r "<< r << " colliding " << is_in_obstacle(loc_coord, test_sphere_radius, DONT_CHECK_WITH_BLOCKING_SPHERES) << " bias " << inside_sampling_bias << endl;
      if((r > inside_sampling_bias && counter == 1) || is_in_obstacle(loc_coord, test_sphere_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)) break;
    }
 return loc_coord;
}


//if beginning position is not specified, program tries to find an initial position by sampling and checking whether the coordinate is inside protein and whether the probe is colliding with it.
//function argument is ratio by which the sampling region is stretched, so it finds a point closer to the centre of the protein
void find_initial_position(double coordinate_stretch){
  double vertX; double vertY; double vertZ;
  double* loc_coord = new double[3];
  int counter = 0;
    while(true){
      vertX = (((double) rand() / (double) RAND_MAX) * (world_size_x * coordinate_stretch)) + (world_size_x_shift + (world_size_x - (world_size_x * coordinate_stretch)) / 2);
      vertY = (((double) rand() / (double) RAND_MAX) * (world_size_y * coordinate_stretch)) + (world_size_y_shift + (world_size_y - (world_size_y * coordinate_stretch)) / 2);
      vertZ = (((double) rand() / (double) RAND_MAX) * (world_size_z * coordinate_stretch)) + (world_size_z_shift + (world_size_y - (world_size_z * coordinate_stretch)) / 2);
      loc_coord[0] = vertX; loc_coord[1] =  vertY; loc_coord[2] = vertZ;
      //print_vector(loc_coord);
      counter++;
      if(!(is_in_obstacle(loc_coord, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)) && is_in_obstacle(loc_coord, test_sphere_radius / 2, DONT_CHECK_WITH_BLOCKING_SPHERES)) break;
      if(counter > 1000){
        cerr << "initial position not found!" << endl;
        exit(0);
      }
    }
  qix = loc_coord[0];
  qiy = loc_coord[1];
  qiz = loc_coord[2];
  delete [] loc_coord;
}

bool termination_check(double* coordinates){
  return !(is_in_obstacle(coordinates, test_sphere_radius, CHECK_WITH_BLOCKING_SPHERES));
}


bool add_node(shared_ptr<vertex> nearest_neighbor, shared_ptr<vertex> new_vertex){
    new_vertex->set_index(tree_index);
   // cout << "add_node: parent index: " << nearest_neighbor->get_index() << endl;
    //cout << "add node: cur_index: " << new_vertex->get_index() << endl;
    new_vertex->set_parent_pointer(nearest_neighbor);
    nearest_neighbor->add_child_pointer(new_vertex, false); //adds index of
    
    local_priority_kdTree_coordinates->add_node(new_vertex);
    add_to_tree(new_vertex->get_location_coordinates(), tree_index, local_priority_kdTree);

   // cout << "add node: kDtree size " << local_priority_kdTree->size << endl;
   // cout << "add node: tree size " << local_priority_kdTree_coordinates->get_size() << endl;     
    tree_index++;
    if(termination_check(new_vertex->get_location_coordinates())){
      backtrack_process_path(new_vertex);
      return false;       
    }
    return true;
}


void try_add_new_point(MPNN::MultiANN<double>* kDtree_passed, double* sampled_coords, Tree* parent_tree, bool global, bool output){
  if(output){
    //cout << "global " << global << " frame " << get_current_frame() << endl;
    //cout << "it's size " << kDtree_passed->size << endl;
  }

  static int counter = 0;
  //cout << "counter " << counter << endl;
 // cout << "kd tree size " << kDtree_passed->size << endl;
  //cout << "tree size " << parent_tree->get_size() << endl;
  int parent_index = find_nearest_vertex(sampled_coords, kDtree_passed);
  shared_ptr<vertex> parent = parent_tree->get_node_pointer(parent_index); //indices mapping is kept 1:1
  shared_ptr<vertex> new_vertex = make_shared<vertex>(sampled_coords, parent,tree_index, probe_radius, get_current_frame(),get_current_frame(), true);
  //todo - redo so it adds stuff to the trees in the step function, interpolation can lead to different scenarios with one or more new nodes
  step_interpolated(new_vertex, parent, MAX_STEP_DISTANCE, min_step, parent_index);
       //function for adding next point to graph points and sending signals to end iterating, use whichever you want                                                                                     //new_vertex is expected to be deleted if step_success_flag is set false, new vertex is added to                                                                                                     //both local and global coordinate representation, and only to local kd tree, the local kd tree is kept because o                                                                                      //the need to track local kd tree when deleting parts of it  
  if(step_success_flag){
    if(get_current_frame() > 1){
      //should be valid always, so added into TEST case
      if(TESTING_ENABLED){
        if(parent->get_last_valid_frame() != get_current_frame()){
          cout << "TEST ERROR: STEP: parent not valid in current frame!" <<endl;
          exit(0);
        }
      }        
    } 
  }

  counter++;
}

//function for adding points to rr tree until break condition is achieved or program runs out of specified max number of iterations 
void find_path(ofstream* stats_filename){
  for(int i = 0; i < iterations; i++){  
    double* sampled_coords = sample(i);
    if(is_local_tree_initialized){
      try_add_new_point(local_priority_kdTree, sampled_coords, local_priority_kdTree_coordinates, false, false);
    }                                                                                   //the need to track local kd tree when deleting parts of it
    delete [] sampled_coords;

    if(is_finished){ //allows to implement another ending condition then just "running out of iterations", currently not implemented
       break;
     }
   }
}  

//for brute force approach
bool path_check(std::map<int, shared_ptr<vertex>>& tested_path, std::map<int, shared_ptr<vertex>>& reference_path){
  double distances_sum = 0;
  double tested_path_length = 0;
  //beginning point is ignored
  int counter = 0; //for skipping first iteration, because you cannot compute lenght withou prev node
  for(std::map<int, shared_ptr<vertex>>::iterator tested_path_iterator = tested_path.begin(); tested_path_iterator != tested_path.end(); tested_path_iterator++){
    if(counter == 0){counter++; continue;}
    double shortest_distance = DBL_MAX; //initiation of shortest distance of tested path's node to reference path on max value
                                        //searching for the closest one in reference path
    for(std::map<int, shared_ptr<vertex>>::iterator reference_path_iterator = reference_path.begin(); reference_path_iterator != reference_path.end(); reference_path_iterator++){
      double cur_distance = compute_metric_eucleidean(tested_path_iterator->second->get_location_coordinates(), reference_path_iterator->second->get_location_coordinates(), 3);
      if(cur_distance < shortest_distance){shortest_distance = cur_distance;}
    }
    distances_sum += shortest_distance;
    double* prev_coordinates = tested_path[tested_path_iterator->second->get_parent_index()]->get_location_coordinates();
    tested_path_length += compute_metric_eucleidean(tested_path_iterator->second->get_location_coordinates(), prev_coordinates, 3);
  }
  //cout << "check duplication "<<distances_sum / tested_path_length << endl;
  return (distances_sum / tested_path_length) < DUPLICATION_CHECK_CONSTANT ? true : false;
}

int is_duplicated(shared_ptr<Path> tested__path){
  if(use_caver_dupcheck){
    int ret = is_tunnel_duplicated(tested__path, paths);
    if(ret != -1){num_of_duplicates++;}
    return ret;
  }
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  map<int, shared_ptr<vertex>> tested_path = tested__path->get_vertices();

  //for debugging---
  static int counter = 0;
  if(counter > MAX_DUPLICATES || paths.size() >= MAX_TUNNELS){
        is_finished = true;
      }
  //----
  int i = 0;    
  for(std::vector<shared_ptr<Path>>::iterator iterator = paths.begin(); iterator != paths.end(); iterator++){
    bool check1 = path_check(tested_path, (*iterator)->get_vertices());
    bool check2 = path_check((*iterator)->get_vertices(), tested_path);
    if(check1 || check2){
      counter++;
  //    cout << "duplicate " << counter << endl;

      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duplication_check_time += duration_cast<microseconds>( t2 - t1 ).count();

      num_of_duplicates++;
      //cout << "TRIAL -- basic dup check 1" << endl << endl; 
      return i;
    }
  i++;
  }
//  cout << "not duplicates " << endl;
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duplication_check_time += duration_cast<microseconds>( t2 - t1 ).count();
  //cout << "TRIAL -- basic dup check 0" << endl << endl; 
  return -1;
} 



double test = 0;
void add_blocking_sphere(double* coordinates, double radius, shared_ptr<Path> path){
  
  shared_ptr<vertex> endpoint = local_priority_kdTree_coordinates->get_node_pointer(path->get_endpoint_index());
  rebuild_blocking_spheres_structure(coordinates, radius); 
 //---delete invalidated tree part 
 while(true){ 
  //parent  
  shared_ptr<vertex> parent = local_priority_kdTree_coordinates->get_node_pointer(local_priority_kdTree_coordinates->get_parent_index(endpoint->get_index())); 
  //check 
  if(!is_in_obstacle(parent->get_location_coordinates(), probe_radius, CHECK_ONLY_BLOCKING_SPHERES) || parent->get_index() == 0){ 
    break; 
  } 
  //new run  
  endpoint = parent; 
 } 
 cut_subtree_in_main_tree(endpoint, false, true);




  //print_vector(coordinates);
  //out << "radius " << radius << endl;
 //cout << "dist " << compute_metric_eucleidean(path->get_beginning_node()->get_location_coordinates(), coordinates) << endl;
 

 //static int counter_post = 0;
 //static int counter_pre = 0;
 //cout << endl << "before " << local_priority_kdTree_coordinates->get_size() << endl;
 //counter_pre += print_colliding(local_priority_kdTree_coordinates->get_vertices());
 //cout << RESETED_TREE_MODE << endl;
 high_resolution_clock::time_point t1 = high_resolution_clock::now(); 
 //if(!RESETED_TREE_MODE) purge_nodes_in_blocking_spheres(path);
 high_resolution_clock::time_point t2 = high_resolution_clock::now();
 test += duration_cast<microseconds>( t2 - t1 ).count(); 

 //cout << "after " << local_priority_kdTree_coordinates->get_size() << endl << endl;
 //counter_post += print_colliding(local_priority_kdTree_coordinates->get_vertices());
 //cout << "pre count " << counter_pre << endl;
 //cout << "post count " << counter_post << endl;
 //cout << "end " << endl << endl;
 //cout << "test " << test << endl;
 //--- 
}

void label_path_nodes(shared_ptr<Path> path, int path_index){
  for(std::map<int, shared_ptr<vertex>>::iterator iterator = path->get_vertices().begin(); iterator != path->get_vertices().end(); iterator++){
    iterator->second->set_in_path(path_index);    
  }
}

void label_used_nodes(shared_ptr<Path> path, int path_index){
 // print_map_indices(path->get_vertices(), "debug");
 // cout << "labeling, frame " << get_current_frame() << endl;
  if(get_current_frame() > 1){
 //   local_priority_kdTree_coordinates->print_nodes();
 //   global_tree_points->print_nodes();
  }
  for(std::map<int, shared_ptr<vertex>>::iterator iterator = path->get_vertices().begin(); iterator != path->get_vertices().end(); iterator++){
  //  cout <<"in main program " << iterator->second->index << " "<< path->ending_index << endl;
    if(get_current_frame() > 1){
   //   cout << "labeling, is path node local? " << iterator->second->is_local() << "path node index " << iterator->second->get_index() << endl;
    }
    iterator->second->is_local() ? local_priority_kdTree_coordinates->get_node_pointer(iterator->second->get_index())->set_in_path(path_index) : global_tree_points->get_node_pointer(iterator->second->get_index())->set_in_path(path_index);
    
  }
}


shared_ptr<vertex> step(shared_ptr<vertex> new_vertex, shared_ptr<vertex> nearest_neighbor_vertex, double goal_distance, 
  double collision_check_step, int parent_index){
 // print_vertex_coordinates(destination);
  static int counter = 0;
  
  //compute differential coordinates
	double x_differential = new_vertex->get_location_coordinates()[0] - nearest_neighbor_vertex->get_location_coordinates()[0]; double kx = signum(x_differential); //program is using differential coordinates instead of absolute ones
	double y_differential = new_vertex->get_location_coordinates()[1] - nearest_neighbor_vertex->get_location_coordinates()[1]; double ky = signum(y_differential); 
	double z_differential = new_vertex->get_location_coordinates()[2] - nearest_neighbor_vertex->get_location_coordinates()[2]; double kz = signum(z_differential); 
    


	double l = compute_metric_eucleidean(new_vertex->get_location_coordinates(),nearest_neighbor_vertex->get_location_coordinates(), 3);        //distance of generated point(source) from nearest atom(destination)
	double l_planar = compute_metric_eucleidean(new_vertex->get_location_coordinates(),nearest_neighbor_vertex->get_location_coordinates(), 2); //distance of generated point from nearest atom in xy projection
	double l_from_source_to_goal = l - goal_distance;                                                               //distance by which the generated point is moved towards the nearest atom
  double alfa = atan(abs(z_differential) / l_planar);                                                             //angle of vertical and horizontal part of source to destiation trajectory
  double beta; x_differential != 0 ? beta = atan(abs(y_differential) / abs(x_differential)) : beta = 0;           //angle of planar parts 

  double new_x; double new_y; double new_z;


  //basic version without any interpolation

  //if the distance is larger than the max step length, the node is moved to the max step 
  
	if(l > goal_distance){
	 //geometric computation of new coordiantes which satisfy the maximal distance of generated point to nearest atom in rr tree
	 new_x = x_differential - kx * cos(alfa) * cos(beta) * l_from_source_to_goal;
	 new_y = y_differential - ky * cos(alfa) * sin(beta) * l_from_source_to_goal;
	 new_z = z_differential - kz * sin(alfa) * l_from_source_to_goal;
	 new_vertex->get_location_coordinates()[0] = new_x + nearest_neighbor_vertex->get_location_coordinates()[0]; new_vertex->get_location_coordinates()[1] = new_y + nearest_neighbor_vertex->get_location_coordinates()[1]; new_vertex->get_location_coordinates()[2] = new_z + nearest_neighbor_vertex->get_location_coordinates()[2];
  }
  if(is_in_obstacle(new_vertex->get_location_coordinates(), probe_radius, CHECK_WITH_BLOCKING_SPHERES)){
    //static int deleted = 1;
    //cout << "del " << deleted << endl;
    //deleted++;
    step_success_flag = false;  
    return NULL; //discard colliding point
  }
  
  
  
  step_success_flag = true;
  return new_vertex;
}

//use simpler way
array<double, 3> move_in_direction(double distance, double alfa, double beta, int kx, int ky, int kz, array<double, 3> old_coordinates){
  double new_x = old_coordinates[0] - kx * cos(alfa) * cos(beta) * distance;
  double new_y = old_coordinates[1] - ky * cos(alfa) * sin(beta) * distance;
  double new_z = old_coordinates[2] - kz * sin(alfa) * distance;
  array<double,3> ret{new_x, new_y, new_z};
  return ret;
}

void interpolate_segment(shared_ptr<vertex> new_vertex, shared_ptr<vertex> nearest_neighbor_vertex, double step, double alfa, double beta, int kx, int ky, int kz){
  std::vector<shared_ptr<vertex>> interpolated_nodes;
  int first_valid = 0;
  int counter = 0;
  bool is_new_first_valid = true;
  double distance = compute_metric_eucleidean(new_vertex->get_location_coordinates(),nearest_neighbor_vertex->get_location_coordinates(), 3);
  array<double, 3> new_coordinates{new_vertex->get_location_coordinates()[0], new_vertex->get_location_coordinates()[1], new_vertex->get_location_coordinates()[2]};
  //cout << "interpolation: distance: before" << distance << endl;
  while(distance > step){
    

    //finds new coordinates and adds them with the referential vertex    
    new_coordinates = move_in_direction(step, alfa, beta, kx, ky, kz, new_coordinates);

    //cout << "interpolation: distance " << distance << endl;
    //cout << "colliding ? " << is_in_obstacle(new_coordinates, probe_radius, CHECK_WITH_BLOCKING_SPHERES) << endl;

    if(!(is_in_obstacle(new_coordinates, probe_radius, CHECK_WITH_BLOCKING_SPHERES))){
      double* coords = copy_array_to_coordinate_pointer(new_coordinates);
      interpolated_nodes.push_back(make_shared<vertex>(coords, 0, probe_radius, get_current_frame(), get_current_frame(), true));
      delete [] coords;
      if(is_new_first_valid){
        first_valid = interpolated_nodes.size() - 1;
        is_new_first_valid = false;
      }  
    } else {
      is_new_first_valid = true;
    }
    counter++;
    distance = compute_metric_eucleidean(new_coordinates, copy_coordinate_pointer_to_array(nearest_neighbor_vertex->get_location_coordinates()));
    if(TESTING_ENABLED){
      //this means that it already should have gone beyond the original point, yet the loop still runs
      if(counter > (MAX_STEP_DISTANCE / step)){
        cout << "TEST ERROR: move_in_direction: interpolation is not working!" << endl;
        exit(0);
      }
    }
  }

  //this means that node closest to parent one is colliding, therefore they are all discarded
  if(is_new_first_valid){
    step_success_flag = false;
    return;
  }

  for(int i = interpolated_nodes.size() - 1; i >= first_valid; i--){
    double epsilon = 0.00000001;
    ////cout << "first valid " << first_valid << endl;
    //cout << "number of nodes " << interpolated_nodes.size() << endl;
    if(i == (interpolated_nodes.size() - 1)){

    //cout << "distance " << compute_metric_eucleidean(nearest_neighbor_vertex->get_location_coordinates(), interpolated_nodes[i]->get_location_coordinates()) << endl;
      if(TESTING_ENABLED){
        //rounding error checked
        if(compute_metric_eucleidean(nearest_neighbor_vertex->get_location_coordinates(), interpolated_nodes[i]->get_location_coordinates()) > (step + 0.001)){
          cout << "TEST ERROR, interpolation: main parent node is in larger distance than step distance!" << endl;
          cout << "distance " << compute_metric_eucleidean(nearest_neighbor_vertex->get_location_coordinates(), interpolated_nodes[i]->get_location_coordinates()) << endl;
          cout << "step " << step << endl;
          exit(0);
        }
      }
      //path is backtracked, trees are reconstructed -> process must stop and begin in new space.
      //this is required because this method tends to use nodes which have virtually the same coordinates as their parents, that is unnecesary and can lead to problems later 
      if(compute_metric_eucleidean(nearest_neighbor_vertex->get_location_coordinates(), interpolated_nodes[i]->get_location_coordinates()) < epsilon){
        i--;
        if(i < first_valid){
          break;
        }
      }
      //cout << "distance " << compute_metric_eucleidean(nearest_neighbor_vertex->get_location_coordinates(), interpolated_nodes[i]->get_location_coordinates()) << endl;
      if(!add_node(nearest_neighbor_vertex,interpolated_nodes[i])){
        break;
      }
    } else {
      if(TESTING_ENABLED){
        if(!double_equals(compute_metric_eucleidean(interpolated_nodes[i + 1]->get_location_coordinates(), interpolated_nodes[i]->get_location_coordinates()), step)){
          cout << "TEST ERROR, interpolation: parent node is in not in step distance!" << endl;
          cout << "distance " << compute_metric_eucleidean(interpolated_nodes[i + 1]->get_location_coordinates(), interpolated_nodes[i]->get_location_coordinates()) << endl;
          cout << "step " << step << endl;
          exit(0);
        }
      }
      if(!add_node(interpolated_nodes[i + 1],interpolated_nodes[i])){
        break;
      }
    }
  }
  step_success_flag = true;
}


//new version
//parameters - interpolation_step, max_step
//rewrite this for fucks sake, it's stupidly complicated for such a simple thing
void step_interpolated(shared_ptr<vertex> new_vertex, shared_ptr<vertex> nearest_neighbor_vertex, double max_step, 
  double interpolation_step, int parent_index){
 // print_vertex_coordinates(destination);
  static int counter = 0;
  
  //compute differential coordinates
  double x_differential = new_vertex->get_location_coordinates()[0] - nearest_neighbor_vertex->get_location_coordinates()[0]; double kx = signum(x_differential); //program is using differential coordinates instead of absolute ones
  double y_differential = new_vertex->get_location_coordinates()[1] - nearest_neighbor_vertex->get_location_coordinates()[1]; double ky = signum(y_differential); 
  double z_differential = new_vertex->get_location_coordinates()[2] - nearest_neighbor_vertex->get_location_coordinates()[2]; double kz = signum(z_differential); 
    


  double l = compute_metric_eucleidean(new_vertex->get_location_coordinates(),nearest_neighbor_vertex->get_location_coordinates(), 3);        //distance of generated point(source) from nearest atom(destination)
  double l_planar = compute_metric_eucleidean(new_vertex->get_location_coordinates(),nearest_neighbor_vertex->get_location_coordinates(), 2); //distance of generated point from nearest atom in xy projection
  double l_from_source_to_goal = l - max_step;                                                               //distance by which the generated point is moved towards the nearest atom
  double alfa = atan(abs(z_differential) / l_planar);                                                             //angle of vertical and horizontal part of source to destiation trajectory
  double beta; x_differential != 0 ? beta = atan(abs(y_differential) / abs(x_differential)) : beta = 0;           //angle of planar parts 

  double new_x; double new_y; double new_z;


  //basic version without any interpolation

  //if the distance is larger than the max step length, the node is moved to the max step 
  //cout << "step_interpolated: original distance " << l << endl;
  if(l > max_step){
   //geometric computation of new coordiantes which satisfy the maximal distance of generated point to nearest atom in rr tree
   new_x = x_differential - kx * cos(alfa) * cos(beta) * l_from_source_to_goal;
   new_y = y_differential - ky * cos(alfa) * sin(beta) * l_from_source_to_goal;
   new_z = z_differential - kz * sin(alfa) * l_from_source_to_goal;
   new_vertex->get_location_coordinates()[0] = new_x + nearest_neighbor_vertex->get_location_coordinates()[0]; new_vertex->get_location_coordinates()[1] = new_y + nearest_neighbor_vertex->get_location_coordinates()[1]; new_vertex->get_location_coordinates()[2] = new_z + nearest_neighbor_vertex->get_location_coordinates()[2];
  }
  if(l > interpolation_step){
    interpolate_segment(new_vertex, nearest_neighbor_vertex, min_step, alfa, beta, kx, ky, kz);
  }
}
  

void add_to_path_count(int path){
  if(!path_found_in_frame[path]){
    paths_count[path] += 1;
    path_found_in_frame[path] = true;
  }
}

void restart(){
  delete local_priority_kdTree;
  delete local_priority_kdTree_coordinates;
  build_tree(LOCAL);
  local_priority_kdTree_coordinates = new Tree();
  double location_coordinates[3] = {qix, qiy, qiz};
  shared_ptr<vertex> start = make_shared<vertex>(location_coordinates, 0, probe_radius, get_current_frame(), get_current_frame(), true);
  local_priority_kdTree_coordinates->add_node(start);
  
  tree_index = 0;  //use as a static variable in tree::add_node;
  add_to_tree(start->get_location_coordinates(), tree_index, local_priority_kdTree);
  tree_index++; 
  is_local_tree_initialized = true;
}

void backtrack_process_path(shared_ptr<vertex> new_vertex ) {   
		found = true;
    //new vertex is the last node added to terminated path
		shared_ptr<Path> path = backtrack(new_vertex);
    if(path == NULL){
      return;
    }

    int last_index = path_optimization(path);
    int is_duplicate = is_duplicated(path);
    bool add_sphere;
    double* location_coordinates = path->get_node_pointer(path->get_endpoint_index())->get_location_coordinates();
    double blocking_sphere_radius = path->get_node_pointer(path->get_endpoint_index())->get_radius();
    //cout << "block " << blocking_sphere_radius << endl;

    if(is_duplicate == -1){
      paths.push_back(path);
      //print_map_indices(path->get_vertices(), "printing path: ");
      label_used_nodes(path, paths.size() - 1);
      duplicate_count.push_back(0);
      paths_count.push_back(1);
      path_found_in_frame.push_back(true);
      add_cluster_path(path, paths.size() - 1);
      add_sphere = true;
    } else{
      add_cluster_path(path, is_duplicate);
      add_to_path_count(is_duplicate);
      duplicate_count[is_duplicate] += 1;   
      add_sphere = true;     
    }

    if(is_duplicate == -1){ 
     add_blocking_sphere(location_coordinates, blocking_sphere_radius, path); 
    } else { 
      if(0.2*duplicate_count[is_duplicate] < maximum_increase){ 
        add_blocking_sphere(location_coordinates, (0.2*duplicate_count[is_duplicate] + 1) * blocking_sphere_radius, path); 
      } 
      else { 
        add_blocking_sphere(location_coordinates, (maximum_increase) * blocking_sphere_radius, path); 
      } 
    } 
 
    if(RESETED_TREE_MODE){
      restart();
    } 
	}


MPNN::MultiANN<double>* new_frame_initalisation(){
  
  //*for clustering
  //marks all paths as not found in new frame, therefore they can be added to counts of the same path in new frame if found again
  for (int i = 0; i < path_found_in_frame.size(); ++i) {
    path_found_in_frame[i] = false;
  }

  run_next_frame(); //switches frame number and loads next coordinate structure

  for (map<int, shared_ptr<vertex>>::iterator iterator = local_priority_kdTree_coordinates->get_vertices().begin(); iterator != local_priority_kdTree_coordinates->get_vertices().end(); iterator++){
    shared_ptr<vertex> cur = iterator->second;
    //use another data structure for inactives to avoid unncessary checks
    if(!cur->is_inactive() && !is_in_obstacle(cur->get_location_coordinates(), probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)){
      iterator->second->set_last_valid_frame(get_current_frame());
    } else {
      iterator->second->set_inactive(true);
    }
    /*
    cout << "index " << iterator->second->get_index() << endl;
    cout << "first frame " << iterator->second->get_first_frame() << endl;
    cout << "last frame " << iterator->second->get_last_valid_frame() << endl;
    cout << "inactive " << iterator->second->is_inactive() << endl;
  */
  }

  //cout << "kd_tree before " << local_priority_kdTree->size << endl;
  rebuild_kd_tree(true, true, true);                  //builds new kd tree without invalid nodes in current frame
  //cout << "kd_tree after " << local_priority_kdTree->size << endl;
}



int find_nearest_vertex(double* loc_coord, MPNN::MultiANN<double>* kdTree){
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  int idx;
  double dann = INFINITY;

  MPNN::ANNpoint query = MPNN::annAllocPt(dimension);
  query[0] = loc_coord[0];
  query[1] = loc_coord[1];
  query[2] = loc_coord[2];
  // find nearest neighbor to the point 'query'. The result is 'nearest' which is index to 'data'
  int ret = kdTree->NearestNeighbor(query,idx,dann);
  MPNN::annDeallocPt(query);

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  nearest_neighbor_search_time += duration_cast<microseconds>( t2 - t1 ).count();
  nearest_vertex_search_count++;

  return ret;
}





int Ball::counter = 0;
int vertex::vertex_counter = 0;


int main(int argc, char *argv[])
{  

  struct stat buf;
  const char* iterations_stats_file_normal = "iterations_stats_normal.txt";
  int exists = stat(iterations_stats_file_normal, &buf);
  if(exists == 0){
    remove(iterations_stats_file_normal);
  }
  std::ofstream num_iterations_stats_normal;
  num_iterations_stats_normal.open(iterations_stats_file_normal); 

  
  load_parameters("config.txt");
  init_protein_struct();
  init_random();
  if(!is_init_position_given()){
    find_initial_position(0.3);
  }
  is_init_position_given() ? cout  << "init position from config file: " << qix << " " << qiy << " " << qiz << endl:cout  << "valid init position found: " << qix << " " << qiy << " " << qiz << endl;
  cout << get_current_frame() << endl;
  cout << "probe radius: " << probe_radius << endl;
  cout << "inside smapling bias " << inside_sampling_bias << endl;
  init_start();
  

  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  for(int i = 0; i < get_frames_count(); i++){
    if(i > 0) {
      new_frame_initalisation();
    }
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    print_frame_borders();
  //---main part
    find_path(&num_iterations_stats_normal);
  //---
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    cout << "time spent on frame: " << get_current_frame() << " : " << duration_cast<microseconds>( t2 - t1 ).count() << endl;
  }

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto main_running_time = duration_cast<microseconds>( t2 - t1 ).count();


  int found_tunnels_count = paths.size();
  if(are_duplicated_pdbs_created){
   //  paths.insert(paths.end(), duplicated_paths.begin(), duplicated_paths.end());
  }

  for(int i = 0; i < paths.size(); i++){
    if(!test_path_noncolliding_static(paths[i])){
      cout << "path colliding before processing!" <<endl;
      create_segfault();
    }
    path_optimization_postprocessing(paths[i]);
    if(!test_path_noncolliding_static(paths[i])){
      cout << "path colliding before processing!" <<endl;
      create_segfault();
    }
  }

  if(found){
    string directory = "./found_tunnels/tunnels";
    if(argc > 1){
      string directory_name(directory + argv[1]);
      directory_name.append("/trajectory");
      write_paths_to_pdbs(directory_name,paths);
    } else {
      string directory_name(directory + "0");
      directory_name.append("/trajectory");
      write_paths_to_pdbs(directory_name,paths);
    }
    if(get_current_frame() > 1){
      string clusters_directory = "./clusters/clusters";  // 0/clusters";
      if(argc > 1){
        clusters_directory.append(argv[1]);
      } else {
        clusters_directory.append("0");
      }
      clusters_directory.append("/cluster");
      write_clusters_to_pdbs(clusters_directory, duplicated_paths);
    }
  }
  //cout << "copies count " << copies_count << endl;
  std::ofstream blocking_file;
  blocking_file.open("blocking_spheres.pdb");
  write_to_pdb(blocking_file, blocking_spheres_in_vertices);
  blocking_file.close();
  
  
  int tree_size_at_finish = local_priority_kdTree_coordinates->get_size();
  
  num_iterations_stats_normal.close();

  stringstream sout;
  sout << "tree size at end: " << tree_size_at_finish << endl;
  sout << "found tunnels count: " << found_tunnels_count << endl;
  sout << "detected dupplicate paths count: " << num_of_duplicates << endl;
  sout << "main part running time: " << main_running_time / 1000000 << endl;
  sout << "time spent searching nearest neighbor: " << nearest_neighbor_search_time / 1000000<< endl;
  sout << "time spent optimizing tunnels: " << tunnel_optimization_time / 1000000<< endl;
  sout << "time spent centring nodes: " << centring_time / 1000000<< endl;
  sout << "time spent collision checking: " << collision_check_time / 1000000<< endl;
  sout << "time spent checking duplications: " << duplication_check_time / 1000000<< endl;
  sout << "time spent purging: " << test / 1000000<< endl;
  sout << endl;
  for(int i = 0; i < paths_count.size(); i++){
    sout << "tunnel "<<  i + 1 << " found in " << ((double)paths_count[i] / (double)get_current_frame())*100 << "% of tested frames" << endl;
  } 

  //create_segfault();
  string data = sout.str();

  cout << data;

  std::ofstream outfile;

  outfile.open("runtime_stats.txt", std::ios_base::app);
  outfile << data;
  outfile << endl;
  

  return 0;
}

