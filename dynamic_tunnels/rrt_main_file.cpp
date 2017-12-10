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
#include "./rapid/RAPID.H"
#include "obstacle.h"
#include "MPNN/include/DNN/multiann.h"
#include "MPNN/include/DNN/ANN.h"
#include "support_methods.h"
#include "path_duplication_check.h"

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
std::vector <shared_ptr<Path>> paths;        //here are the fpund paths kept, all nodes in them are copied, therefore there isn't problem with deleting
                                  //their nodes everywhere. However, that doesn't protect you, if you want to access the node in main structure!
                                  //All nodes have related indices in main structure, however there is no guarantee that the ones in main structure
                                  //weren't deleted!

std::vector <int> duplicate_count; //keeps track of path count in "clusters"

//--- used to keep track for tunnel validity in frames
vector <int> paths_count;
vector <bool> path_found_in_frame;
//-----

Tree* local_priority_kdTree_coordinates; //exists only because dumbass MPNN doesn't support delete, kepts coordinates of local kd Tree
Tree* global_tree_points;                //to:do
double maximum_increase = 1.5; //aximum increase in blocking sphere size due to many duplicates, in percents (1.5 = 150 percent bigger sphere)
bool step_success_flag;        //used to determine succesful add to local tree, if fails global is tryied or if global doesn't exist, iteration fail
bool is_local_tree_initialized;
bool LOCALIZED_MODE = false; 
//kept for debugging purposes---
bool are_duplicated_pdbs_created = false;
std::map <int, shared_ptr<vertex>> blocking_spheres_in_vertices;
std::vector <shared_ptr<Path>> duplicated_paths;
//vector <vertex*> visualise_centring;
//---

MPNN::MultiANN<double> *global_kdTree;
MPNN::MultiANN<double> *local_priority_kdTree;

//-----variables for finding program statistics
auto nearest_neighbor_search_time = 0;
auto tunnel_optimization_time = 0;
auto centring_time = 0;
auto collision_check_time = 0;
auto duplication_check_time = 0;
int nearest_vertex_search_count = 0;
int num_of_duplicates = 0;
//-----

const int dimension = 3;                     //k-d tree dimension

//------program parameters  //consider using txt file
const double DUPLICATION_CHECK_CONSTANT = 0.6; //ratio of inter-tunnel distance to tunnel length, if lower, tunnel is considered a duplicate of existing tunnel
const double BLOCKING_SPHERE_RADIUS_COEFFICIENT = 1.2;
//******these are used for debugging purposes, set them to some absurd number for normal run
const int MAX_DUPLICATES = 10000;               //num of duplicates before run terminates //
const int MAX_TUNNELS = 10000;                   //maximum number of tunnels before run terminates
//******
//------

//--- identificators of global/local kd_tree for easier reading
const int NONE = 0;
const int GLOBAL = 1; 
const int LOCAL = 2;
const int BOTH = 3;
//---

int tree_index = 0;


//obvious
shared_ptr<vertex> load_parent_vertex_from_tree(shared_ptr<vertex> cur_vertex){
  return cur_vertex->is_parent_local() ? local_priority_kdTree_coordinates->get_node_pointer(cur_vertex->get_parent_index()) : global_tree_points->get_node_pointer(cur_vertex->get_parent_index()); 
}

//frame is valid if both tested node and it's child are valid in it and and the tested node's is same or smaller "is previous"
//when locating parent node, we can't "travel into the future", the node should be in the same or older frame than the actual. If it is valid in more than one,obviously the most recent is used.
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

//sets node into already found frame, deletes all children other than child in path, returns pointer to node added to path
//deletes other children pointers so the path has only the path child, 
shared_ptr<vertex> copy_node_tree_to_path(shared_ptr<Path> path, shared_ptr<vertex> tree_node, shared_ptr<vertex> path_child_node, int cur_frame){


  shared_ptr<vertex> path_cur_node = tree_node->copy_without_structure_pointers();
  path_cur_node->set_frame_index(cur_frame);
  path_cur_node->add_child_pointer(path_child_node, false);
  path_child_node->set_parent_pointer(path_cur_node);
  path->add_node(path_cur_node);



  return path_cur_node;
}

//creates a copy of vertices involved in path from found endpoint to initial point, the path must be copied in case that it's part is deleted
shared_ptr<Path> backtrack(shared_ptr<vertex> endpoint){ 
  shared_ptr<Path> path = make_shared<Path>(-1, endpoint->get_index(), get_current_frame());


  shared_ptr<vertex> cur = load_parent_vertex_from_tree(endpoint);
  cout << "Tree: backtrack: parent index " << endpoint->get_parent_index() << endl;

  if(cur->get_last_valid_frame() != endpoint->get_last_valid_frame()) {cerr << "error in backtrack, endpoint and it's parents last frames doesn't match! " << endl; exit(0);} //these are more for testing, shoudn't happen
  if(endpoint->get_last_valid_frame() != get_current_frame()) {cerr << "error in backtrack, endpoint isn't in current frame! " << endl; exit(0);}
  if(endpoint->get_children_count() > 0) {cerr << "error in backtrack, endpoint has children! " << endl; exit(0);}

  shared_ptr<vertex> endpoint_copied = endpoint->copy_without_structure_pointers(); //
  endpoint_copied->set_frame_index(get_current_frame());
  path->add_node(endpoint_copied);

  shared_ptr<vertex> child_pointer = copy_node_tree_to_path(path, cur, endpoint_copied, get_current_frame());
  

  while(cur->get_parent_index() != -1){ //cur is previous node from the tree
    cur = load_parent_vertex_from_tree(cur);
    cout << "backtrack: parent index " << cur->get_parent_index() << endl;

    int nearest_frame = find_valid_frame(cur, child_pointer);
    
    if(nearest_frame > child_pointer->get_frame_index()){
      cout << "backtrack_frame cur " << nearest_frame << endl;
      cout << "backtrack_frame prev " << child_pointer->get_frame_index() << endl;
      cout << "backtrack fail, time travel not allowed!" << endl;
    }
    if(nearest_frame == -1){            //there is no previous frame, therefore you would have time travel to past to go through this path
      cerr << "found null frame in path " << endl; //SHOULDN'T HAPPEN!
      sleep(10);
      return NULL;
    }
    child_pointer = copy_node_tree_to_path(path, cur, child_pointer, nearest_frame);
  }
  path->set_beginning_index(cur->get_index()); 
  
  return path;
}

//initializes k-d tree, puts initial vertex with parent index -1, symbolizing it hasn't got a parent, into tree map
void init_start(){
  build_tree(LOCAL);
  local_priority_kdTree_coordinates = new Tree();
  double location_coordinates[3] = {qix, qiy, qiz};
  shared_ptr<vertex> start = make_shared<vertex>(location_coordinates, 0, probe_radius, get_current_frame(), true);
  local_priority_kdTree_coordinates->add_node(start);
  
  tree_index = 0;  //use as a static variable in tree::add_node;
  add_to_tree(start->get_location_coordinates(), tree_index, local_priority_kdTree);
  tree_index++; 
  found = false;
  is_local_tree_initialized = true;
}

//collision checking without checking the blocking spheres structure is used in functions not involving path search (path optimization and chcecking whether new node is inside the protein)
const int CHECK_ONLY_BLOCKING_SPHERES = 2;
const int CHECK_WITH_BLOCKING_SPHERES = 1;
const int DONT_CHECK_WITH_BLOCKING_SPHERES = 0;
bool is_in_obstacle(double* loc_coord, double radius, int check_with_blocking_spheres){
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  static int col_counter = 0;
  ozcollide::Sphere sphere;
  sphere.center.x = loc_coord[0];
  sphere.center.y = loc_coord[1];
  sphere.center.z = loc_coord[2];
  sphere.radius = radius;

  bool returned;

  high_resolution_clock::time_point t2;

  switch ( check_with_blocking_spheres ) {
    case CHECK_ONLY_BLOCKING_SPHERES:
    if(num_blocking_spheres() >= 1){
      AABBTreeSphere* blocking_spheres_tree = get_blocking_spheres_tree();

      t2 = high_resolution_clock::now();
      collision_check_time += duration_cast<microseconds>( t2 - t1 ).count();

      return blocking_spheres_tree->isCollideWithSphere(sphere);   
    } else {

      t2 = high_resolution_clock::now();
      collision_check_time += duration_cast<microseconds>( t2 - t1 ).count();

      return false;}
    break;
    case DONT_CHECK_WITH_BLOCKING_SPHERES:
      returned = protein_tree->isCollideWithSphere(sphere);
      t2 = high_resolution_clock::now();
      collision_check_time += duration_cast<microseconds>( t2 - t1 ).count();
      return returned;
    break;
    case CHECK_WITH_BLOCKING_SPHERES:
     bool has_collided_with_protein_structure = protein_tree->isCollideWithSphere(sphere);
     bool has_collided_with_blocking_spheres;
     if(num_blocking_spheres() >= 1){
        AABBTreeSphere* blocking_spheres_tree = get_blocking_spheres_tree();
        has_collided_with_blocking_spheres = blocking_spheres_tree->isCollideWithSphere(sphere);   
      } else {has_collided_with_blocking_spheres = false;}

     t2 = high_resolution_clock::now();
     collision_check_time += duration_cast<microseconds>( t2 - t1 ).count();

     return (has_collided_with_protein_structure || has_collided_with_blocking_spheres);
    break;
  }
}

//lazy hack
bool is_in_obstacle(array<double, 3> loc_coord, double radius, int check_with_blocking_spheres){
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  static int col_counter = 0;
  ozcollide::Sphere sphere;
  sphere.center.x = loc_coord[0];
  sphere.center.y = loc_coord[1];
  sphere.center.z = loc_coord[2];
  sphere.radius = radius;

  bool returned;

  high_resolution_clock::time_point t2;

  switch ( check_with_blocking_spheres ) {
    case CHECK_ONLY_BLOCKING_SPHERES:
    if(num_blocking_spheres() >= 1){
      AABBTreeSphere* blocking_spheres_tree = get_blocking_spheres_tree();

      t2 = high_resolution_clock::now();
      collision_check_time += duration_cast<microseconds>( t2 - t1 ).count();

      return blocking_spheres_tree->isCollideWithSphere(sphere);   
    } else {

      t2 = high_resolution_clock::now();
      collision_check_time += duration_cast<microseconds>( t2 - t1 ).count();

      return false;}
    break;
    case DONT_CHECK_WITH_BLOCKING_SPHERES:
      returned = protein_tree->isCollideWithSphere(sphere);
      t2 = high_resolution_clock::now();
      collision_check_time += duration_cast<microseconds>( t2 - t1 ).count();
      return returned;
    break;
    case CHECK_WITH_BLOCKING_SPHERES:
     bool has_collided_with_protein_structure = protein_tree->isCollideWithSphere(sphere);
     bool has_collided_with_blocking_spheres;
     if(num_blocking_spheres() >= 1){
        AABBTreeSphere* blocking_spheres_tree = get_blocking_spheres_tree();
        has_collided_with_blocking_spheres = blocking_spheres_tree->isCollideWithSphere(sphere);   
      } else {has_collided_with_blocking_spheres = false;}

     t2 = high_resolution_clock::now();
     collision_check_time += duration_cast<microseconds>( t2 - t1 ).count();

     return (has_collided_with_protein_structure || has_collided_with_blocking_spheres);
    break;
  }
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
    cout << "add_node: parent index: " << nearest_neighbor->get_index() << endl;
    cout << "add node: cur_index: " << new_vertex->get_index() << endl;
    new_vertex->set_parent_pointer(nearest_neighbor);
    nearest_neighbor->add_child_pointer(new_vertex, false); //adds index of
    
    local_priority_kdTree_coordinates->add_node(new_vertex);
    add_to_tree(new_vertex->get_location_coordinates(), tree_index, local_priority_kdTree);

    cout << "add node: kDtree size " << local_priority_kdTree->size << endl;
    cout << "add node: tree size " << local_priority_kdTree_coordinates->get_size() << endl;     
    tree_index++;
    if(termination_check(new_vertex->get_location_coordinates())){
      backtrack_process_path(new_vertex);
      return false;       
    }
    return true;
}

//kDtree_passed - kDtree in which the nearest point is looked for, global ??
void try_add_new_point(MPNN::MultiANN<double>* kDtree_passed, double* sampled_coords, Tree* parent_tree, bool global, bool output){
  if(output){
    cout << "global " << global << " frame " << get_current_frame() << endl;
    cout << "it's size " << kDtree_passed->size << endl;
  }

  static int counter = 0;
  cout << "counter " << counter << endl;
  cout << "kd tree size " << kDtree_passed->size << endl;
  cout << "tree size " << parent_tree->get_size() << endl;
  int parent_index = find_nearest_vertex(sampled_coords, kDtree_passed);
  shared_ptr<vertex> parent = parent_tree->get_node_pointer(parent_index); //indices mapping is kept 1:1
  shared_ptr<vertex> new_vertex = make_shared<vertex>(sampled_coords, parent,tree_index, probe_radius, get_current_frame(), true);
  //todo - redo so it adds stuff to the trees in the step function, interpolation can lead to different scenarios with one or more new nodes
  step_interpolated(new_vertex, parent, MAX_STEP_DISTANCE, min_step, parent_index);
       //function for adding next point to graph points and sending signals to end iterating, use whichever you want                                                                                     //new_vertex is expected to be deleted if step_success_flag is set false, new vertex is added to                                                                                                     //both local and global coordinate representation, and only to local kd tree, the local kd tree is kept because o                                                                                      //the need to track local kd tree when deleting parts of it  
  if(step_success_flag){
    if(get_current_frame() > 1){
      //should be valid always, so added into TEST case
      if(TESTING_ENABLED){
        if(!parent->is_valid_in_frame(get_current_frame())){
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
      }
    if(LOCALIZED_MODE && get_current_frame() > 1 && (!step_success_flag || !is_local_tree_initialized)){
      try_add_new_point(global_kdTree, sampled_coords, global_tree_points, true, false);
      if(step_success_flag){
        is_local_tree_initialized = true;
      }
    }                                                                                     //the need to track local kd tree when deleting parts of it
    delete [] sampled_coords;

    if(is_finished){ //allows to implement another ending condition then just "running out of iterations", currently not implemented
       break;
     }
   }
}  

//function which returns approximation of maximal radius of sphere fitting into the tunnel, expects non-colliding sphere 
double approximate_radius(double* coordinates, double default_radius, double precision, double default_approximation_step){
  double approximation_step = default_approximation_step;
  double cur_radius = default_radius + approximation_step;
  while(approximation_step > precision && cur_radius < test_sphere_radius + 1){
    bool hasCollided = is_in_obstacle(coordinates, cur_radius, DONT_CHECK_WITH_BLOCKING_SPHERES);
    if(!hasCollided){
      cur_radius += approximation_step;
    }
    else{
      cur_radius -= approximation_step;
      approximation_step = approximation_step / 2;
      cur_radius += approximation_step;
    }
  }
  return cur_radius;
}


/*
Function which returns node approximately centered around the tunnel path. It can be used in two modes, either it centers node without respect to neighboring nodes, which leads to jagged path
but is more accurate or it takes into account two neighboring nodes and moving them in the same direction as the main node aimed for centering, which results in smoother path but more approximate
centering. Function returns new node without deleting the original, so the user must delete the original one if that is desired. Neighboring nodes' coordinates are changed only.
Function firslty finds plane perpendicular to supplied direction vector, using two perpendicular vectors to it and then performs "gradient descent" trying to move the node in the plane in attempt
to find the maximum possible radius to fit inside the tunnel. 
*/
double* center_node(shared_ptr<vertex> centered_node, shared_ptr<vertex> previous_node, shared_ptr<vertex> next_node, double* coordinate_vector, double increment, bool are_neighbors_centered, bool rewrite_node, bool return_node_info){
  high_resolution_clock::time_point t1 = high_resolution_clock::now();


  double normalised_length = 1;
  //what must be deleted perpendicular vector, node_coordinates, cur_perpendicular_vector, cur_second_perpendicular_vector, shift_vector, shifted_coordinates
  //first perpendicular vector, found by transposition of coordinates
  double perpendicular_vector_2D1[3]  = {-coordinate_vector[1], coordinate_vector[0], 0};
  double perpendicular_vector_2D2[3]  = {0, coordinate_vector[2], -coordinate_vector[1]};
  double* perpendicular_vector = add_vectors(perpendicular_vector_2D1, perpendicular_vector_2D2, ADDITION);
  
  //second base perpendicular vector, found using cross product
  double second_perpendicular_vector[3]  = {coordinate_vector[1]*perpendicular_vector[2] - coordinate_vector[2]*perpendicular_vector[1],
  coordinate_vector[2]*perpendicular_vector[0] - coordinate_vector[0]*perpendicular_vector[2],
  coordinate_vector[0]*perpendicular_vector[1] - coordinate_vector[1]*perpendicular_vector[0]};

  //normalization of perpendicular (plane base) vectors to specified length
  normalize_vector(perpendicular_vector, normalised_length);
  normalize_vector(second_perpendicular_vector, normalised_length);
  
  //coordinates of centered node, it must be copied because the function is meant to keep the original node unchanged. Neighboring are declared and intialized to use later in case neighbor centering
  //is allowed 
  double* node_coordinates = centered_node->copy_coordinates();
  double* prev_node_shifted_coordinates;
  double* next_node_shifted_coordinates;
  if(are_neighbors_centered){
    prev_node_shifted_coordinates = previous_node->copy_coordinates();
    if(next_node != NULL){
      next_node_shifted_coordinates = next_node->copy_coordinates(); //happens for last node, that one should be centered, unlike the root node
    }
      
  } 
  
  //maximum radius in default position, used in first iteration to find larger one 
  double cur_radius = approximate_radius(node_coordinates, centered_node->get_radius(), 0.01, 0.1);
  //for comparison
  double check_radius = cur_radius; 
  double default_radius = centered_node->get_radius(); //reference radius, important because radius adjusting function expects non-colliding molecule, if
  //radius used to store previous maximum value to check the one found in current iteration against 
  double prev_radius = default_radius;

  /*  //sanity check
  cout << "--dot products--" << endl;
  cout << dot_product(coordinate_vector, perpendicular_vector) << endl;
  cout << dot_product(coordinate_vector, second_perpendicular_vector) << endl;
  cout << dot_product(second_perpendicular_vector, perpendicular_vector) << endl;
  cout << "--dot products out--" << endl;
  */

  //gradient descent, in a plane bounded by two perpendicular plane vectors, because we don't want to find places outside the nodes postion in tunnel
  //debug **
  /*
  static int counter = 0;      //**
  counter++;                   //**
  static int vis_counter = 0;  //** 
  */

  //parameters of algorithm, min_increment is the minimum shift of node at which the algortihm still runs and number_of_trials is amount of different translations
  //the algorithm tries in one iteration
  double min_increment = 0.125;
  int number_of_trials = 10;

  while((increment > min_increment) && cur_radius < test_sphere_radius){
    //variable to store current best increase in maximum radius which will be found by next iteration of algorithm
    double max_increase = 0;
    //runs number_of_trials interations with translations of fixed length and variable directions in translation plane from current node position with attempt to find new position with largest radius
    for(int i = 0; i < number_of_trials; i++){
      //angle according to which the node will be shifted in this attempt
      double angle = (i * 2 * M_PI) / (number_of_trials - 1);                  //deterministic descent
      //double angle = ((double) rand() / (double) RAND_MAX) * 2 * M_PI;       //stochastic descent

      //decomposition of generated angle into base vectors
      double first_vector_part = increment * normalised_length * sin(angle);   
      double second_vector_part = increment * normalised_length * cos(angle);
      
      //base vectors are copied, so the original ones wouldn't be lost by shifting by current angle and multiplied by base vector components of current tested direction
      double* cur_perpendicular_vector = copy_vector(perpendicular_vector); double* cur_second_perpendicular_vector = copy_vector(second_perpendicular_vector);
      multiply_vector(cur_perpendicular_vector, first_vector_part); multiply_vector(cur_second_perpendicular_vector, second_vector_part);
      //finally the shifted vectors are added to point to location in circle with increment radius and summed with the node's coordinates
      double* shift_vector = add_vectors(cur_perpendicular_vector, cur_second_perpendicular_vector, ADDITION);      
      double* shifted_coordinates = add_vectors(shift_vector, node_coordinates, ADDITION);

      //... and maximum radius of ball in that placement is found
      cur_radius = approximate_radius(shifted_coordinates, default_radius, 0.01, 0.1);
     
      //maximum radius in this position is checked subtracted from maximum radius found in previous iteration and checked against the maximum increase found in current iteration, which starts at zero,
      //therefore the new maximum radius must be greater than the previous one
      if((cur_radius - prev_radius) > max_increase && !is_in_obstacle(shifted_coordinates, cur_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)) {
        

        //--- can be buggy, check later
        if(are_neighbors_centered){
          //new coordinates shifted in the same way as main node, this is to prevent excessive deviation and resulting jaggedness of tunnel
          double* neighbors_shift_vector = copy_vector(shift_vector); multiply_vector(neighbors_shift_vector, 0.5);
         // cout << neighbors_shift_vector[0] << endl;
         // cout << shift_vector[0] << endl;
          double* cur_prev_node_shifted_coordinates = add_vectors(neighbors_shift_vector, prev_node_shifted_coordinates, ADDITION);
          double* cur_next_node_shifted_coordinates = add_vectors(neighbors_shift_vector, next_node_shifted_coordinates, ADDITION);
          delete [] neighbors_shift_vector;

          if((next_node == NULL && !is_in_obstacle(prev_node_shifted_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES))||(next_node != NULL && !is_in_obstacle(prev_node_shifted_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES) && !is_in_obstacle(next_node_shifted_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES))){
             //test was succesfull, node can be moved
             //new radius is assigned
             max_increase = cur_radius - prev_radius;
             prev_radius = cur_radius;
             delete [] node_coordinates;
             node_coordinates = copy_vector(shifted_coordinates);
             //if replacement is possible, the old translated nodes' coordinates are replaced with new one
             delete [] prev_node_shifted_coordinates; delete [] next_node_shifted_coordinates; 
             prev_node_shifted_coordinates = copy_vector(cur_prev_node_shifted_coordinates);
             next_node_shifted_coordinates = copy_vector(cur_next_node_shifted_coordinates);
             delete [] cur_prev_node_shifted_coordinates; delete [] cur_next_node_shifted_coordinates;
          } 
          else {
             delete [] cur_prev_node_shifted_coordinates; delete [] cur_next_node_shifted_coordinates;
          }
        }
        

        //---can be buggy. check later---end
        else {
          max_increase = cur_radius - prev_radius;
          prev_radius = cur_radius;
          delete [] node_coordinates;
          node_coordinates = copy_vector(shifted_coordinates); 
       // cout << cur_radius << endl;
        }
  
      }
      delete [] cur_perpendicular_vector; delete [] cur_second_perpendicular_vector; delete [] shift_vector; delete [] shifted_coordinates;
   }
   //cout << endl;
   increment = increment / 1.5;

  }
  /*
  if(counter % 20 == 0){vis_counter++; cout << "centring field size " << visualise_centring.size()<<endl;} //**
  */
  double ret_radius;
  cur_radius > check_radius ? ret_radius = cur_radius : ret_radius = check_radius;
  delete [] perpendicular_vector;
  if(rewrite_node){
  centered_node->set_coordinates(node_coordinates);
  centered_node->set_radius(cur_radius);
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  centring_time += duration_cast<microseconds>( t2 - t1 ).count();

  if(are_neighbors_centered){
    delete [] prev_node_shifted_coordinates;
    if(next_node != NULL){
      delete [] next_node_shifted_coordinates; //happens for last node, that one should be centered, unlike the root node
    }    
  } 


  if(!return_node_info){ //this is because we don't always want to allocate new memory
    delete [] node_coordinates;
    return NULL;
  } else {
    double* ret = new double[4];
    ret[0] = node_coordinates[0]; ret[1] = node_coordinates[1]; ret[2] = node_coordinates[2]; ret[3] = cur_radius;
    delete [] node_coordinates; 
    return ret;   
  }
  //std::cout<< cur_radius << endl;
}


void smoothing(shared_ptr<Path> path){
  int i = 0;
  int next_atom_index = 2;
  int beginning_index = path->get_beginning_index(); 
  //cout << path.size() << endl;

  std::map<int, shared_ptr<vertex>>& path_vertices = path->get_vertices();
  shared_ptr<vertex> cur = path_vertices[beginning_index];
  if(path_vertices.size() < 3){
    return;
  }

  

  int cur_index = cur->get_index();
  int child_index = path->get_child_index(cur_index);
  int second_child_index = path->get_child_index(child_index);
  while(path->get_children_count(second_child_index) > 0){
    cur_index = cur->get_index();
    child_index = path->get_child_index(cur_index);
    second_child_index = path->get_child_index(child_index);    //if crap like this is needed more in the future, consider modifying path structure to index nodes within path for random access inside


    double* direction_vector = add_vectors(path_vertices[second_child_index]->get_location_coordinates(), path_vertices[cur_index]->get_location_coordinates(), SUBTRACTION);
    multiply_vector(direction_vector, 0.5);
    double* shifted_coordinates_vector = add_vectors(direction_vector, path_vertices[cur_index]->get_location_coordinates(), ADDITION);
    bool is_not_passable = is_in_obstacle(shifted_coordinates_vector, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES);
    if(!is_not_passable){
      path->erase_node_with_reconnecting(child_index);
    } 

    cur = path_vertices[second_child_index];
    //path->print_nodes();
    cur_index = cur->get_index();
    //cout << "cur_index " << cur_index;
    child_index = path->get_child_index(cur_index);
    //cout << "child_index " << child_index;
    if(path_vertices[child_index]->get_children_count() == 0){ //can happen because of deleting passable elements
      delete [] shifted_coordinates_vector;
      delete [] direction_vector;
      break;
    }
    second_child_index = path->get_child_index(child_index);
    //cout << "scond child_index " << second_child_index;

    delete [] shifted_coordinates_vector;
    delete [] direction_vector;
  }
}

//Later I should research how to use generics here instead of two functions
void rebuild_kd_tree(bool is_initialized, int status, bool copy_stuff, bool only_valid_in_curframe){
  cout << "rebuild status " << status <<endl;
  if(is_initialized){
    if(status == GLOBAL || status == BOTH) {delete global_kdTree;}
    if(status == LOCAL || status == BOTH) {delete local_priority_kdTree;}
  }
  if(status == GLOBAL || status == BOTH)build_tree(GLOBAL);
  if(status == LOCAL || status == BOTH)build_tree(LOCAL);
  
  if(copy_stuff && (status == GLOBAL || status == BOTH)){
    std::map<int, shared_ptr<vertex>> graph_points = global_tree_points->get_vertices();
    cout << "rebuild glob " << graph_points.size() << endl;
    for (std::map<int, shared_ptr<vertex>>::iterator iterator = graph_points.begin(); iterator != graph_points.end(); iterator++){
      add_to_tree(iterator->second->get_location_coordinates(), iterator->second->get_index(), global_kdTree);
    }
  }
  if(copy_stuff && (status == LOCAL || status == BOTH)){

    std::map<int, shared_ptr<vertex>> graph_points_loc = local_priority_kdTree_coordinates->get_vertices();
    cout << "rebuild loc " << graph_points_loc.size() << endl;
    for (std::map<int, shared_ptr<vertex>>::iterator iterator = graph_points_loc.begin(); iterator != graph_points_loc.end(); iterator++){
      if(!only_valid_in_curframe){
        add_to_tree(iterator->second->get_location_coordinates(), iterator->second->get_index(), local_priority_kdTree);
      }
      else{
        if(iterator->second->get_last_valid_frame() == get_current_frame()){
          add_to_tree(iterator->second->get_location_coordinates(), iterator->second->get_index(), local_priority_kdTree);
        }
      } 
    }
    cout << "comparison: kd tree " << local_priority_kdTree->size << " data struct " << graph_points_loc.size() << endl;  
  }
}




//recursive function deleting all children of particular node, remember paths nodes are copied and have their respective children at ending position of chidren_vertices(.. i should
//delete the previous.. to_do .. ), but indices are matching with the main graph_parts structure.. deleting paths and subtrees should be thought of as separate functions. Prev index is imporatant
//for recuresive part of the run for deleting child references in parent nodes, therefore we don't care about it in original call
void cut_subtree_in_main_tree_recursion(shared_ptr<vertex> cur, bool has_debugging_output, bool is_local){
  //cout << "recuresion, children count " <<cur->get_children_count() << endl;
  while(cur->get_children_count() != 0){
    /*
    if(has_debugging_output){
      cout << "--printing path--" << endl;
      for(std::map<int, shared_ptr<vertex>>::iterator iterator = graph_points.begin(); iterator != graph_points.end(); iterator++){
        cout << iterator->second->get_index() << " "; 
      }
      cout << endl <<"--printing path--end--" << endl;
      cout << "index "<< cur->children_indices[cur->children_indices.size() -1] << endl;
      cout << "bef "<< cur->children_indices.size() << endl; 
    }
    */
    shared_ptr<vertex> next = cur->get_child_pointer(cur->get_children_count() -1);
    cut_subtree_in_main_tree_recursion(next, has_debugging_output, next->is_local());
    if(has_debugging_output){
      cout << "aft "<< cur->get_children_count() << endl; 
    }
    cur->delete_child(cur->get_children_count() -1); //last one remaining
  }
  if(has_debugging_output){cout << "deleting " <<cur->get_index()<<endl;}
  cur->is_local() ?  local_priority_kdTree_coordinates->delete_node(cur->get_index()) : global_tree_points->delete_node(cur->get_index());

}

//deletes subtree beginning at index and deletes reference to subtree root from it's parent node 
void cut_subtree_in_main_tree(shared_ptr<vertex> first_delet_node, bool has_debugging_output){

  shared_ptr<vertex> node_in_tree;
  shared_ptr<vertex> parent_node_in_tree;

  //cout << "cut: nodes: " << endl;
 // local_priority_kdTree_coordinates->print_nodes();
  if(get_current_frame() > 1){
  //  global_tree_points->print_nodes();
  }
  //cout << "is first local? " << first_delet_node->is_local() << endl;
  first_delet_node->is_local() ? node_in_tree = local_priority_kdTree_coordinates->get_node_pointer(first_delet_node->get_index()) : node_in_tree = global_tree_points->get_node_pointer(first_delet_node->get_index());
  parent_node_in_tree = node_in_tree->get_parent_pointer();
  bool found = false;
  //cout << "path_index " << first_delet_node->get_index() << " path parent index "<< parent_path->get_index()<<" tree_index " << node_in_tree->get_index() << endl;
  for(int i = 0; i < parent_node_in_tree->get_children_count(); i++){
    if(parent_node_in_tree->get_child_index(i) == node_in_tree->get_index()){
      parent_node_in_tree->delete_child(i);
      found = true;
      break;
    }
  }
  if(found == false){
    cout << "size of local tree " << local_priority_kdTree_coordinates->get_size() << endl;
    cout << "size of global tree " << global_tree_points->get_size() << endl;
    cout << "path node index " << first_delet_node->get_index() << endl;
    cout << "parent tree node children count " << endl;
    cerr << "fatal error - in cutting subtree, root of subtree doesn't have parent" << endl;
    exit(0);
  }

  int global_tree_size;
  global_tree_points != NULL ? global_tree_size = global_tree_points->get_size() : 0;
  int local_tree_size = local_priority_kdTree_coordinates->get_size();
  int kdTree_rebuild_index = 0; //four value logic, need to find which tree(s) to rebuild, 0 - none, 1 - global, 2 - local, 3 - both
  int prev = local_priority_kdTree_coordinates->get_size();
  cut_subtree_in_main_tree_recursion(node_in_tree, has_debugging_output, first_delet_node->is_local());

  if(global_tree_points != NULL){
    if(global_tree_points->get_size() != global_tree_size)kdTree_rebuild_index++; 
  }
  if(local_priority_kdTree_coordinates->get_size() != local_tree_size)kdTree_rebuild_index += 2;
  rebuild_kd_tree(true, kdTree_rebuild_index, true, true);
  cout << "tree rebuilt" << endl;
}


int copies_count = 0;
void center_tunnel(shared_ptr<Path> path){
  std::map<int, shared_ptr<vertex>>::iterator iterator;

  int counter = 0;
  int index = -1;
 
  std::map<int, shared_ptr<vertex>>& path_vertices = path->get_vertices();
  cout << "center_tunnel: path_vertices size " << path_vertices.size() << endl;
  for(iterator = path_vertices.begin(); iterator != path_vertices.end(); iterator++){
    //start node is not centered
    if(counter == 0){counter++;index = path_vertices[0]->get_child_index(); continue;}
    
    /*
    With reused tree, that already centered node is inside some previous trajectory. Therefore instead of costy centering procedure it can just be copied from that path.
    exists_in_path() method returns -1 if the node is new. else it returns index into path in paths vector.
    Important - node must be valide in the same frame as current node, because molecular positions are different in every frame!
    Duplicated paths are deleted, therefore we must be sure to not index into them, this is achieved by calling the "labeling" function, which points the nodes in main data structures
    into valid paths, when the path is added into main paths vector 
    */
    shared_ptr<vertex> cur_node = path_vertices[index];
    bool is_local = cur_node->is_local();

    int exists_in_path = -1;
    if(is_local){
      exists_in_path = local_priority_kdTree_coordinates->get_node_pointer(cur_node->get_index())->is_in_path();
    } else {
      exists_in_path = global_tree_points->get_node_pointer(cur_node->get_index())->is_in_path();
    }

    if(exists_in_path != -1){
      shared_ptr<vertex> related_node = paths[exists_in_path]->get_node_pointer(cur_node->get_index());
      cur_node->set_coordinates(related_node->get_location_coordinates()); //copy identical nodes coordinates to current;
      cur_node->set_radius(related_node->get_radius());
      copies_count++;

      if(cur_node->get_children_count() > 0){
        index = cur_node->get_child_index();
      }
      continue;
    }


    shared_ptr<vertex> prev_node = path_vertices[cur_node->get_parent_index()];
    cout << "center tunnel: counter " << counter << " cur_node " << cur_node->get_index() << endl;
    shared_ptr<vertex> child_node = path_vertices[path->get_child_index(cur_node->get_index())];

    double* direction_vector = add_vectors(cur_node->get_location_coordinates(), prev_node->get_location_coordinates(), SUBTRACTION);
    
    center_node(cur_node, prev_node, child_node, direction_vector, 1, true, true, false);
  
    delete [] direction_vector;
    
    /* deleting subtrees is important because of path blocking. Blocking sphere needs to be placed in the center of the tunnel at it's endpoint
    and the tree after this endpoint would create would create valid growth into the empty region, which obviously sucks. It's
    caused by the fact that the part of the tunnel growing after original endpoint ends after centering in empty space.
    both local and global coordinates of the tree after endpoint needs to be declared invalid now and deleted from kd trees,
    and the kd trees are rebuilded. In first frame, it's possible to only care about local structures, which speeds up static detection.
    */



    if(path_vertices[index]->get_radius() > test_sphere_radius){
      int child_index = path->get_child_index(index);    
      if(child_index == -1){
          cerr << "your start position already lies in empty space!" << endl;
          exit(0);      
      }
      cut_subtree_in_main_tree(path->get_node_pointer(child_index), false); 
      path->erase_from_index(child_index); 
      break;}; //to_do
    if(cur_node->get_children_count() > 0){
      index = cur_node->get_child_index(0);
    }
    counter++;
  }
  
  path->set_endpoint_index(index); 
}




void cut_tunnel(shared_ptr<Path> path){
  std::map<int, shared_ptr<vertex>>& path_vertices = path->get_vertices();
  std::vector<int> vertex_indices_to_be_deleted;

  int endpoint_index = path->get_endpoint_index();
  int endpoint_parent_index = path->get_parent_index(endpoint_index);
  //print_map_indices(path->get_vertices(), "cut tunnel");
  while((endpoint_index != -1 && endpoint_parent_index != -1) && path->get_node_pointer(endpoint_index)->is_in_path() == -1){   //check also with global
    //---
    /*
    path->print_nodes();
    cout << endl << "ep " << endpoint_index;
    cout << endl << "epp " << endpoint_parent_index << endl; 
    //---
    */
    //cout << "cut tunnel indices " << endpoint_index << " " << endpoint_parent_index << endl;
    //cout << "cut tunnel " << path_vertices[endpoint_index]->index << " " << path_vertices[endpoint_parent_index]->index << " " << local_priority_kdTree_coordinates[endpoint_index]->is_in_path()<< endl;
    double* direction_vector = add_vectors(path_vertices[endpoint_index]->get_location_coordinates(), path_vertices[endpoint_parent_index]->get_location_coordinates(), SUBTRACTION);
    double* node_information = center_node(path_vertices[endpoint_index], NULL, NULL, direction_vector, test_sphere_radius / 1.5, false, false, true);
    if(node_information [3] > test_sphere_radius){      
      vertex_indices_to_be_deleted.push_back(endpoint_index);
    } 

    delete [] direction_vector;
    delete [] node_information;

    endpoint_index = path->get_parent_index(endpoint_index);
    endpoint_parent_index = path->get_parent_index(endpoint_index);

  }

  int deleted_index = path->get_child_index(vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1]);
  shared_ptr<vertex> cur = path->get_node_pointer(vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1]);

  int prev = path->get_size();
  cut_subtree_in_main_tree(path->get_node_pointer(deleted_index),false); //from element before the last
  path->erase_from_index(deleted_index);
  //cout << "path post " << prev - path->get_size() << endl;
  //exit(0);
  
}

//function checking whether tested path duplicates part of reference path, complete duplication check must be two way

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
      cout << "TRIAL -- basic dup check 1" << endl << endl; 
      return i;
    }
  i++;
  }
//  cout << "not duplicates " << endl;
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duplication_check_time += duration_cast<microseconds>( t2 - t1 ).count();
  cout << "TRIAL -- basic dup check 0" << endl << endl; 
  return -1;
} 

//returns index of last element in cut path.
int path_optimization(shared_ptr<Path> path){
  //path = smoothing(path);
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
 
  smoothing(path);
  center_tunnel(path);
  cut_tunnel(path);
  

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  tunnel_optimization_time += duration_cast<microseconds>( t2 - t1 ).count();
  
  return 1;
  //cut_tunnel(path);
  //path = center_tunnel(path);
}

void add_blocking_sphere(double* coordinates, double radius){
  static int blocking_spheres_counter = 0;
 shared_ptr<vertex> blocking_sphere = make_shared<vertex>(coordinates, 0, radius, get_current_frame(), true); //don't care about index
 blocking_spheres_in_vertices.insert(std::pair<int,shared_ptr<vertex>>(blocking_spheres_counter,blocking_sphere));
 blocking_spheres_counter++;
 rebuild_blocking_spheres_structure(coordinates, radius);
 

}

void label_used_nodes(shared_ptr<Path> path, int path_index){
 // print_map_indices(path->get_vertices(), "debug");
  cout << "labeling, frame " << get_current_frame() << endl;
  if(get_current_frame() > 1){
 //   local_priority_kdTree_coordinates->print_nodes();
 //   global_tree_points->print_nodes();
  }
  for(std::map<int, shared_ptr<vertex>>::iterator iterator = path->get_vertices().begin(); iterator != path->get_vertices().end(); iterator++){
  //  cout <<"in main program " << iterator->second->index << " "<< path->ending_index << endl;
    if(get_current_frame() > 1){
      cout << "labeling, is path node local? " << iterator->second->is_local() << "path node index " << iterator->second->get_index() << endl;
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
  cout << "interpolation: distance: before" << distance << endl;
  while(distance > step){
    

    //finds new coordinates and adds them with the referential vertex    
    new_coordinates = move_in_direction(step, alfa, beta, kx, ky, kz, new_coordinates);

    cout << "interpolation: distance " << distance << endl;
    cout << "colliding ? " << is_in_obstacle(new_coordinates, probe_radius, CHECK_WITH_BLOCKING_SPHERES) << endl;

    if(!(is_in_obstacle(new_coordinates, probe_radius, CHECK_WITH_BLOCKING_SPHERES))){
      interpolated_nodes.push_back(make_shared<vertex>(copy_array_to_coordinate_pointer(new_coordinates), 0, probe_radius, get_current_frame(), true));
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

    cout << "first valid " << first_valid << endl;
    cout << "number of nodes " << interpolated_nodes.size() << endl;

    if(i == (interpolated_nodes.size() - 1)){
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
  cout << "step_interpolated: original distance " << l << endl;
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
  //all shit is put into local kdtree first, and on frame change the global tree is reconstructed from rrt
  /*
  add_to_tree(new_vertex->location_coordinates, tree_index, local_priority_kdTree); //all shit is put in local tree first, and on frame change
  
  tree_index++;
  
	if(termination_check(new_vertex->location_coordinates))
  */

void add_to_path_count(int path){
  if(!path_found_in_frame[path]){
    paths_count[path] += 1;
    path_found_in_frame[path] = true;
  }
}
void backtrack_process_path(shared_ptr<vertex> new_vertex ) {    
		found = true;
    //new vertex is the last node added to terminated path
		shared_ptr<Path> path = backtrack(new_vertex);
    if(path == NULL){
      return;
    }

    path_optimization(path);
    int is_duplicate = is_duplicated(path);

    bool add_sphere;
       
    //---finds strongly centered version of endpoint index, it's radius and coordinates
    /*
    std::map<int, shared_ptr<vertex>>& path_vertices = path->get_vertices();
    double* direction_vector = add_vectors(path_vertices[path->ending_index]->location_coordinates, path_vertices[path->get_parent_index(path->ending_index)]->location_coordinates, SUBTRACTION);
    double* centered_node_info = center_node(path->get_node_pointer(path->ending_index), NULL, NULL, direction_vector, test_sphere_radius, false, false, true); delete [] direction_vector;
    double* location_coordinates = new double[3]; location_coordinates[0] = centered_node_info[0]; location_coordinates[1] = centered_node_info[1]; location_coordinates[2] = centered_node_info[2];
    double blocking_sphere_radius = centered_node_info[3];
    delete [] centered_node_info;
    */
    //---

    double* location_coordinates = path->get_node_pointer(path->get_endpoint_index())->get_location_coordinates();
    double blocking_sphere_radius = path->get_node_pointer(path->get_endpoint_index())->get_radius();

    if(is_duplicate == -1){
      paths.push_back(path);
      //print_map_indices(path->get_vertices(), "printing path: ");
      label_used_nodes(path, paths.size() - 1);
      duplicate_count.push_back(0);
      paths_count.push_back(1);
      path_found_in_frame.push_back(true);
      add_sphere = true;
    } else{
      if(are_duplicated_pdbs_created){
        duplicated_paths.push_back(path);
        add_to_path_count(is_duplicate);
        duplicate_count[is_duplicate] += 1;
        cout << paths_count[is_duplicate] << endl;
        //exit(0);
        add_sphere = true;
      }
      else {
        add_to_path_count(is_duplicate);
        cout << paths_count[is_duplicate] << endl;
        duplicate_count[is_duplicate] += 1;
        add_sphere = true;
      }
    }
    if(is_duplicate == -1){
     add_blocking_sphere(location_coordinates, blocking_sphere_radius);
    } else {
      if(0.2*duplicate_count[is_duplicate] < maximum_increase){
        add_blocking_sphere(location_coordinates, (0.2*duplicate_count[is_duplicate] + 1) * blocking_sphere_radius);
      }
      else {
        add_blocking_sphere(location_coordinates, (maximum_increase) * blocking_sphere_radius);
      }
    }   
	}

void print_vertex_coordinates(shared_ptr<vertex> vert){
	std::cout <<"printing vertex location coordinates: x: " << vert->get_location_coordinates()[0] << " y: " << vert->get_location_coordinates()[1] << " z : " << vert->get_location_coordinates()[2] << endl 
	<< "parent_index: " << vert->get_parent_index() << endl;
}

void copy_to_global(){
  for(std::map<int, shared_ptr<vertex>>::iterator iterator = local_priority_kdTree_coordinates->get_vertices().begin(); iterator != local_priority_kdTree_coordinates->get_vertices().end(); iterator++){
    global_tree_points->add_node(iterator->second);
    iterator->second->make_global();
  }
  for(int i = 0; i < paths.size(); i++){
    for(std::map<int, shared_ptr<vertex>>::iterator iterator = paths[i]->get_vertices().begin(); iterator != paths[i]->get_vertices().end(); iterator++){
      iterator->second->make_global();
    }
  }
}

MPNN::MultiANN<double>* new_frame_initalisation(){
  if(get_current_frame() == 1 && LOCALIZED_MODE){
    global_tree_points = new Tree();
  }
  for (int i = 0; i < path_found_in_frame.size(); ++i) {
    path_found_in_frame[i] = false;
  }
  run_next_frame(); //switches frame number and loads next coordinate structure
  if(LOCALIZED_MODE){
    copy_to_global();
    local_priority_kdTree_coordinates->reset(); //deletes storage for nodes in previous local tree
  }
  //checks which vertices from the old frame are not colliding in current frame, marks them as valid ones

  if(LOCALIZED_MODE){
    for (map<int, shared_ptr<vertex>>::iterator iterator = global_tree_points->get_vertices().begin(); iterator != global_tree_points->get_vertices().end(); iterator++){
      if(!is_in_obstacle(iterator->second->get_location_coordinates(), probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)){
        iterator->second->add_valid_frame(get_current_frame());
      }
    }
  } else {
    for (map<int, shared_ptr<vertex>>::iterator iterator = local_priority_kdTree_coordinates->get_vertices().begin(); iterator != local_priority_kdTree_coordinates->get_vertices().end(); iterator++){
      if(!is_in_obstacle(iterator->second->get_location_coordinates(), probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)){
        iterator->second->add_valid_frame(get_current_frame());
      }
    }
  }
  int i = 1;
  //checks endpoints of paths, if they are valid in current frame. If succeeds, the path is automatically valid in this frame too, cause it
  //always possible to connect it to previous' frame's one
  delete_blocking_spheres();
  for (vector<shared_ptr<Path>>::iterator iterator = paths.begin(); iterator != paths.end(); iterator++){
    if(!is_in_obstacle((*iterator)->get_node_coordinates((*iterator)->get_endpoint_index()), test_sphere_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)){
      (*iterator)->add_valid_frame(get_current_frame());
      add_blocking_sphere((*iterator)->get_node_coordinates((*iterator)->get_endpoint_index()), (*iterator)->get_node_pointer((*iterator)->get_endpoint_index())->get_radius()); //to:do save in vector instead and rebuild it only once!, lazy and tired now
      add_to_path_count(i - 1);
    } else {
      cout << "path " << i << " is invalid in frame " << get_current_frame() << endl;
    }
    i++;
  }
  //first frame doesn't have priority tree, therefore it is not needed to delete it in second frame
  //first bool labels. whether old kd_tree exists, and hence is delet, second which kd_tree is rebuilded, and third whether it's corresponding sturcuture is copied into it
  //locals are not rebuilded, all points in them are copied into new globals 
  if(get_current_frame() == 2 && LOCALIZED_MODE){
    rebuild_kd_tree(false, GLOBAL, true, false);
  } else if(LOCALIZED_MODE) {
    rebuild_kd_tree(true, GLOBAL, true, false);
  }

  if(LOCALIZED_MODE){
    rebuild_kd_tree(true, LOCAL, false, false);                //deletes last frame's priority tree
    is_local_tree_initialized = false;
  }
  if(!LOCALIZED_MODE){
    cout << "kd_tree before " << local_priority_kdTree->size << endl;
    rebuild_kd_tree(true, LOCAL, true, true);                  //builds new kd tree without invalid nodes in current frame
    cout << "kd_tree after " << local_priority_kdTree->size << endl;

  }

     //init priority tree for current frame
}



void build_tree(int is_global){
  
    int *topology = new int[dimension];
    MPNN::ANNcoord *scale = new MPNN::ANNcoord[dimension];
    for(int i=0;i<dimension;i++) {
      scale[i] = 1;
      topology[i] = 1;
 
    }
    if(is_global == GLOBAL){
      global_kdTree = new MPNN::MultiANN<double>(dimension,1,topology,(MPNN::ANNpoint)scale);
    } else {
      local_priority_kdTree = new MPNN::MultiANN<double>(dimension,1,topology,(MPNN::ANNpoint)scale);
    }
    delete [] topology;
    delete [] scale;
}

int add_to_tree(double* loc_coord, int tree_index, MPNN::MultiANN<double>* kdTree){
    //std::cout << "print " <<endl;
    MPNN::ANNpoint annpt = MPNN::annAllocPt(dimension);  
    annpt[0] = loc_coord[0];
    annpt[1] = loc_coord[1];
    annpt[2] = loc_coord[2];
    //print_vector(loc_coord);
    //std::cout << annpt[0] << " " << annpt[1] << " " << annpt[2] << endl;
    kdTree->AddPoint(annpt,tree_index); 
    MPNN::annDeallocPt(annpt);
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



void write_to_pdb(ofstream& file, std::map<int, shared_ptr<vertex>>& atoms){
  std::map<int, shared_ptr<vertex>>::iterator iterator;
  for(iterator = atoms.begin(); iterator != atoms.end(); iterator++){
      file << setprecision(3) << fixed;
      file << "ATOM";

      file.width(7); file << iterator->second->get_frame_index();
      float c1 = 1;
      file << "  N   ILE E  16";
      file.width(12); file <<  iterator->second->get_location_coordinates()[0]; file.width(8); file <<  iterator->second->get_location_coordinates()[1]; file.width(8); file <<  iterator->second->get_location_coordinates()[2]; file.width(2); file.precision(2); file << "  "; file.width(4); file << c1;
      file.width(1); file << " "; file.width(5); file << iterator->second->get_radius(); file.width(12); file << "N"<<endl;      
  }
}


void write_paths_to_pdbs(std::string filename, std::vector<shared_ptr<Path>>& paths){

  int counter = 0;
  std::vector<shared_ptr<Path>>::iterator iterator;
  int i = 0;
  for(iterator = paths.begin(); iterator != paths.end(); iterator++){
    std::ofstream file;
    std::string filename_appended(filename + to_string(i + 1));
    filename_appended.append(".pdb");
    file.open(filename_appended); 
    std::cout << i << " " << filename_appended << endl;
    write_to_pdb(file, (*iterator)->get_vertices());
    file.close();
    i++;
 }
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
    paths.insert(paths.end(), duplicated_paths.begin(), duplicated_paths.end());
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
  }
  cout << "copies count " << copies_count << endl;
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
  sout << endl;
  for(int i = 0; i < paths_count.size(); i++){
    sout << "tunnel "<<  i + 1 << " found in " << ((double)paths_count[i] / (double)get_current_frame())*100 << "% of tested frames" << endl;
  } 

  string data = sout.str();

  cout << data;

  std::ofstream outfile;

  outfile.open("runtime_stats.txt", std::ios_base::app);
  outfile << data;
  outfile << endl;
  
 
  return 0;
}

