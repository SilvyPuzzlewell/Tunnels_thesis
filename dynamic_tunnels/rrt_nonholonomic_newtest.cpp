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
int find_nearest_vertex(double* loc_coord);
void write_paths_to_pdb(const char* filename, vector<vertex*>& path);
double* collision_function(double x1, double y1, double x2, double y2, double min_step_col);
int add_to_tree(double* loc_coord, int tree_index, MPNN::MultiANN<double>* kdTree);
void build_tree();
void print_vertex_coordinates(vertex* vert);
std::map<int, vertex*>& step(std::map <int, vertex*>& graph_points, vertex* generated, vertex* parent, double goal_distance, double collision_check_step, int parent_index);

bool found = false;
bool is_finished = false;
std::vector <Path*> paths;
std::vector <int> duplicate_count;
double maximum_increase = 1.5; //aximum increase in blocking sphere size due to many duplicates, in percents (1.5 = 150 percent bigger sphere)
//kept for debugging purposes---
bool are_duplicated_pdbs_created = false;
std::map <int, vertex*> blocking_spheres_in_vertices;
std::vector <Path*> duplicated_paths;
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
const int INSIDE_SAMPLING_BIAS = 3;          //minimum ratio of samples inside the protein in form n:1
const double DUPLICATION_CHECK_CONSTANT = 0.6; //ratio of inter-tunnel distance to tunnel length, if lower, tunnel is considered a duplicate of existing tunnel
const double BLOCKING_SPHERE_RADIUS_COEFFICIENT = 1.2;
//******these are used for debugging purposes, set them to some absurd number for normal run
const int MAX_DUPLICATES = 10000;               //num of duplicates before run terminates //
const int MAX_TUNNELS = 10000;                   //maximum number of tunnels before run terminates
//******
//------


int tree_index = 0;

void clear_graph_points(map<int, vertex*>&  graph_points){
  std::map<int, vertex*>::iterator iterator;
  for(iterator = graph_points.begin(); iterator != graph_points.end(); iterator++){
    delete iterator->second;
  }
  graph_points.clear(); 
}

//creates a copy of vertices involved in path from found endpoint to initial point
Path* backtrack(map <int, vertex*>& graph_points, int path_end_index){ 
  vertex* cur = graph_points[path_end_index];
  Path* path = new Path(-1, path_end_index, get_current_frame());
  //path.push_back(cur->copy());

  //int child_index = -1;  //marks the node as childless, as last in the path

  vertex* path_node1 = cur->copy();
  path_node1->frame_index = get_current_frame();
  path_node1->children_indices.clear();
  //path_node->children_indices.push_back(child_index);
  path->add_node(cur->index, path_node1);
  cur = path_node1;


  int i = 1;
  while(cur->parent_index != -1){
    //std::cout << "bt " << cur->parent_index << " " << graph_points.size() << endl;
    int child_index = cur->index;
    int child_frame = cur->frame_index;
    cur = graph_points[cur->parent_index];
    // cout << "debug " << child_frame << endl;
    //search for
    int nearest_frame = cur->find_previous_frame(child_frame);
    if(nearest_frame == -1){            //there is no previous frame, therefore you would have time travel to past to go through this path
      delete path;
      return NULL;
    }
    vertex* path_node = cur->copy();
    path_node->children_indices.clear();
    path_node->children_indices.push_back(child_index);
    path_node->frame_index = nearest_frame;
   // path.push_back(cur->copy());
    path->add_node(cur->index, path_node);
    cur = path_node;
    i++;
  }
  path->beginning_index = cur->index; 
  
  return path;
}

//initializes k-d tree, puts initial vertex with parent index -1, symbolizing it hasn't got a parent, into tree map
std::map<int, vertex*>& init_start(std::map <int, vertex*>& graph_points){
  build_tree();
  double location_coordinates[3] = {qix, qiy, qiz};
  int parent_index = -1;
  vertex* start = new vertex(location_coordinates, parent_index, 0, probe_radius, get_current_frame());
  graph_points.insert(std::pair<int,vertex*>(0, start));
  
  tree_index = 0;
  add_to_tree(graph_points[0]->location_coordinates, tree_index, global_kdTree);
  tree_index++; 
  found = false;
  return graph_points;
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
  

//function for generating next sample 
double* sample(int iterations_counter){
  double vertX; double vertY; double vertZ;
  double* loc_coord = new double[3];
    while(true){
      vertX = (((double) rand() / (double) RAND_MAX) * (world_size_x)) + world_size_x_shift;
      vertY = (((double) rand() / (double) RAND_MAX) * (world_size_y)) + world_size_y_shift;
      vertZ = (((double) rand() / (double) RAND_MAX) * (world_size_z)) + world_size_z_shift;
      loc_coord[0] = vertX; loc_coord[1] =  vertY; loc_coord[2] = vertZ;
      if(iterations_counter % INSIDE_SAMPLING_BIAS == 0 || is_in_obstacle(loc_coord, test_sphere_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)) break;
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
//function for adding points to rr tree until break condition is achieved or program runs out of specified max number of iterations 
std::map<int, vertex*>& find_path(std::map<int, vertex*>& graph_points, ofstream* stats_filename){
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  int iterations_counter = 0;

  for(int i = 0; i < iterations; i++){
    iterations_counter++;

    double* sampled_coords = sample(iterations_counter);    
    vertex* new_vertex = new vertex(sampled_coords, parent_index, tree_index, probe_radius, get_current_frame()); //generated point, same index as for accessing the kd_tree is used


    if(get_current_frame() > 1){ 
      int parent_index = find_nearest_vertex(sampled_coords, local_priority_kdTree);
      vertex* parent = graph_points[parent_index];
      graph_points = step(graph_points, new_vertex, parent, min_step, 1, parent_index, local_priority_kdTree); //function for adding next point to graph points and sending signals to end iterating, use whichever you want
    }
    if(step_success_flag == false || get_current_frame == 1){

    }
    delete [] sampled_coords;

    if(is_finished){
       high_resolution_clock::time_point t2 = high_resolution_clock::now();  //code for computing run statistics
       auto duration = duration_cast<microseconds>( t2 - t1 ).count();      
       *stats_filename << " " <<iterations_counter << " " << duration << "  " << graph_points.size() << endl;
       return graph_points;

     }  
  } 

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    *stats_filename << "  " <<iterations_counter << "  " << duration << "  " << graph_points.size() <<  endl;
    return graph_points;
  
}

int signum(double value){
  return value < 0 ? -1 : 1;
}

//definition of function for checking whether goal was found
bool termination_check(double* coordinates){
  return !(is_in_obstacle(coordinates, test_sphere_radius, CHECK_WITH_BLOCKING_SPHERES));
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
double* center_node(vertex* centered_node, vertex* previous_node, vertex* next_node, double* coordinate_vector, double increment, bool are_neighbors_centered, bool rewrite_node, bool return_node_info){
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
  double cur_radius = approximate_radius(node_coordinates, centered_node->radius, 0.01, 0.1);
  //for comparison
  double check_radius = cur_radius; 
  double default_radius = centered_node->radius; //reference radius, important because radius adjusting function expects non-colliding molecule, if
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
      //sanity checks--- **
      /*
      if(counter % 20 == 0 && vis_counter < 50){
        vertex* vert = new vertex(shifted_coordinates, 1, cur_radius);
        visualise_centring.push_back(vert);
      }
      */  
      //---- **

      //cout << "final length " << vector_length(shift_vector) << endl;
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
             prev_node_shifted_coordinates = cur_prev_node_shifted_coordinates;
             next_node_shifted_coordinates = cur_next_node_shifted_coordinates;
          } 
          else {
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
  centered_node->radius = cur_radius;
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  centring_time += duration_cast<microseconds>( t2 - t1 ).count();
  if(!return_node_info){ //this is because we don't always want to allocate new memory
    delete [] node_coordinates;
    return NULL;
  } else {
    double* ret = new double[4];
    ret[0] = node_coordinates[0]; ret[1] = node_coordinates[1]; ret[2] = node_coordinates[2]; ret[3] = cur_radius; 
    return ret;   
  }
  //std::cout<< cur_radius << endl;
}


void smoothing(Path* path){
  int i = 0;
  int next_atom_index = 2;
  int beginning_index = path->beginning_index; 
  //cout << path.size() << endl;

  std::map<int, vertex*>& path_vertices = path->get_vertices();
  vertex* cur = path_vertices[beginning_index];
  int cur_index = cur->index;
  int child_index = path->get_child_index(cur_index);
  int second_child_index = path->get_child_index(child_index);

  //while(path_vertices[second_child_index]->children_indices.size() > 0){
  while(path->get_children_count(second_child_index) > 0){
    cur_index = cur->index;
    child_index = path->get_child_index(cur_index);
    second_child_index = path->get_child_index(child_index);    //if crap like this is needed more in the future, consider modifying path structure to index nodes within path for random access inside


    double* direction_vector = add_vectors(path_vertices[second_child_index]->location_coordinates, path_vertices[cur_index]->location_coordinates, SUBTRACTION);
    multiply_vector(direction_vector, 0.5);
    double* shifted_coordinates_vector = add_vectors(direction_vector, path_vertices[cur_index]->location_coordinates, ADDITION);
    bool is_not_passable = is_in_obstacle(shifted_coordinates_vector, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES);
    if(!is_not_passable){
      path->erase_node_with_reconnecting(child_index);
    } 

    cur = path_vertices[second_child_index];
    //path->print_nodes();
    cur_index = cur->index;
    //cout << "cur_index " << cur_index;
    child_index = path->get_child_index(cur_index);
    //cout << "child_index " << child_index;
    if(path_vertices[child_index]->children_indices.size() == 0){ //can happen because of deleting passable elements
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



//recursive function deleting all children of particular node, remember paths nodes are copied and have their respective children at ending position of chidren_vertices(.. i should
//delete the previous.. to_do .. ), but indices are matching with the main graph_parts structure.. deleting paths and subtrees should be thought of as separate functions. Prev index is imporatant
//for recuresive part of the run for deleting child references in parent nodes, therefore we don't care about it in original call



void rebuild_kd_tree(std::map<int, vertex*> graph_points){
  delete global_kdTree;
  build_tree();
  
  for (std::map<int, vertex*>::iterator iterator = graph_points.begin(); iterator != graph_points.end(); iterator++){
    add_to_tree(iterator->second->location_coordinates, iterator->second->index, global_kdTree);
  }
}
void cut_subtree_in_main_tree_recursion(int index, std::map<int ,vertex*>& graph_points, bool has_debugging_output){
  vertex* cur = graph_points[index];
  while(cur->children_indices.size() != 0){
    //cout << cur->children_indices[cur->children_indices.size() -1] << endl;
    if(has_debugging_output){
      cout << "--printing path--" << endl;
      for(std::map<int, vertex*>::iterator iterator = graph_points.begin(); iterator != graph_points.end(); iterator++){
        cout << iterator->second->index << " "; 
      }
      cout << endl <<"--printing path--end--" << endl;
      cout << "index "<< cur->children_indices[cur->children_indices.size() -1] << endl;
      cout << "bef "<< cur->children_indices.size() << endl; 
    }
    cut_subtree_in_main_tree_recursion(cur->children_indices[cur->children_indices.size() -1], graph_points, has_debugging_output);
    if(has_debugging_output){
      cout << "aft "<< cur->children_indices.size() << endl; 
    }
    cur->children_indices.pop_back();
  }
  if(has_debugging_output){cout << "deleting " <<cur->index<<endl;}
  vertex* deletedptr = graph_points[cur->index];
  graph_points.erase(cur->index);
  delete deletedptr;

}

//deletes subtree beginning at index and deletes reference to subtree root from it's parent node 
void cut_subtree_in_main_tree(int index, std::map<int ,vertex*>& graph_points, bool has_debugging_output){
  int prev_index = graph_points[index]->parent_index;
  vertex* parent = graph_points[prev_index];
  bool found = false;
  for(int i = 0; i < parent->children_indices.size(); i++){
    if(parent->children_indices[i] == index){
      parent->children_indices.erase(parent->children_indices.begin() + i);
      found = true;
      break;
    }
  }
  if(found == false){cerr << "fatal error - in cutting subtree, root of subtree doesn't have parent" << endl; exit(0);}
  cut_subtree_in_main_tree_recursion(index, graph_points, has_debugging_output);
  rebuild_kd_tree(graph_points);
}


//since path has only one child index and empty chidren indices vector denotes endpoint, function for cutting up paths is simpler than trees.
void cut_path(int index, std::map<int ,vertex*>& path){
  vertex* cur = path[index];
  path[cur->parent_index]->children_indices.clear(); //destroy parent reference to first deleted node, thus making parent the path endpoint
  int next_index = -1;
  while(cur->children_indices.size() > 0){
    next_index = cur->children_indices[0];

    vertex* deletedptr = path[cur->index];
    path.erase(cur->index);
    delete deletedptr;

    cur = path[next_index];
  }
  if(next_index == -1){
    cout << "one-atom path!" << endl;
    return;
  }
  vertex* deletedptr = path[cur->index];
  path.erase(cur->index);
  delete deletedptr;
}


int center_tunnel_last_index_flag;
void center_tunnel(Path* path, std::map<int ,vertex*>& graph_points){
  std::map<int, vertex*>::iterator iterator;

  int counter = 0;
  int index = -1;
 /* cout << "path indices " << endl;
  for(std::map<int, vertex*>::iterator iterator = path.begin(); iterator != path.end(); iterator++){
    cout << iterator->second->index << endl;
  }
  cout << "path indices " << endl << endl; */
  std::map<int, vertex*>& path_vertices = path->get_vertices();
  for(iterator = path_vertices.begin(); iterator != path_vertices.end(); iterator++){
    if(counter == 0){counter++;index = path_vertices[0]->children_indices.back(); continue;}

    vertex* cur_node = path_vertices[index];
 //   cout << "parent index " << cur_node->parent_index << endl;
    vertex* prev_node = path_vertices[cur_node->parent_index];
    vertex* child_node = path_vertices[path->get_child_index(cur_node->index)];
 //   cout << "node centering_called" << endl;
    double* direction_vector = add_vectors(cur_node->location_coordinates, prev_node->location_coordinates, SUBTRACTION);
    center_node(cur_node, prev_node, child_node, direction_vector, 1, true, true, false);
    delete [] direction_vector;
    
    
    if(path_vertices[index]->radius > test_sphere_radius){
      int child_index = path->get_child_index(index);
      cut_subtree_in_main_tree(child_index, graph_points, false); 
      cut_path(child_index, path_vertices); 
      break;}; //to_do
    if(cur_node->children_indices.size() > 0){
      index = cur_node->children_indices.back();
    }
  }
  
  path->ending_index = index; 
  center_tunnel_last_index_flag = index; //variable storing last index of center tunnel procedure
}




void cut_tunnel(Path* path, std::map<int, vertex*>& graph_points){
  std::map<int, vertex*>& path_vertices = path->get_vertices();
  std::vector<int> vertex_indices_to_be_deleted;

  int endpoint_index = path->ending_index;
  int endpoint_parent_index = path->get_parent_index(endpoint_index);
  int last_deleted_index = -1; //for deleting subtree in main tree

  while(endpoint_parent_index != -1){
    //---
    /*
    path->print_nodes();
    cout << endl << "ep " << endpoint_index;
    cout << endl << "epp " << endpoint_parent_index << endl; 
    //---
    */
    double* direction_vector = add_vectors(path_vertices[endpoint_index]->location_coordinates, path_vertices[endpoint_parent_index]->location_coordinates, SUBTRACTION);
    double* node_information = center_node(path_vertices[endpoint_index], NULL, NULL, direction_vector, test_sphere_radius, false, false, true);
    if(node_information [3] > test_sphere_radius){

      vertex_indices_to_be_deleted.push_back(endpoint_index);
      last_deleted_index = endpoint_index;

    } else {
      delete [] direction_vector;
      delete [] node_information ;
      break;
    }

    delete [] direction_vector;
    delete [] node_information;

    endpoint_index = path->get_parent_index(endpoint_index);
    endpoint_parent_index = path->get_parent_index(endpoint_index);

  }
  for(int i = 0; i < vertex_indices_to_be_deleted.size() - 1; i++){ //last added element is not deleted, because we want to block the tunnel with it
    path->erase_node_with_reconnecting(vertex_indices_to_be_deleted[i]);
  } 
  if(vertex_indices_to_be_deleted.size() > 1){ //cut if something was deleted from path
    cut_subtree_in_main_tree(vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 2], graph_points, false); //from element before the last
  }
}


int return_min(int i1, int i2){
  return i1 < i2 ? i1 : i2;
}

//function checking whether tested path duplicates part of reference path, complete duplication check must be two way

bool path_check(std::map<int, vertex*>& tested_path, std::map<int, vertex*>& reference_path){
  double distances_sum = 0;
  double tested_path_length = 0;
  //beginning point is ignored
  int counter = 0; //for skipping first iteration, because you cannot compute lenght withou prev node
  for(std::map<int, vertex*>::iterator tested_path_iterator = tested_path.begin(); tested_path_iterator != tested_path.end(); tested_path_iterator++){
    if(counter == 0){counter++; continue;}
    double shortest_distance = DBL_MAX; //initiation of shortest distance of tested path's node to reference path on max value
                                        //searching for the closest one in reference path
    for(std::map<int, vertex*>::iterator reference_path_iterator = reference_path.begin(); reference_path_iterator != reference_path.end(); reference_path_iterator++){
      double cur_distance = compute_metric_eucleidean(tested_path_iterator->second->location_coordinates, reference_path_iterator->second->location_coordinates, 3);
      if(cur_distance < shortest_distance){shortest_distance = cur_distance;}
    }
    distances_sum += shortest_distance;
    double* prev_coordinates = tested_path[tested_path_iterator->second->parent_index]->location_coordinates;
    tested_path_length += compute_metric_eucleidean(tested_path_iterator->second->location_coordinates, prev_coordinates, 3);
  }
  //cout << "check duplication "<<distances_sum / tested_path_length << endl;
  return (distances_sum / tested_path_length) < DUPLICATION_CHECK_CONSTANT ? true : false;
}

int is_duplicated(Path* tested__path){
   cout << "TRIAL -- caver dup check " << is_tunnel_duplicated(tested__path, paths) << endl;
   high_resolution_clock::time_point t1 = high_resolution_clock::now();
   map<int, vertex*> tested_path = tested__path->get_vertices();

  //for debugging---
  static int counter = 0;
  if(counter > MAX_DUPLICATES || paths.size() >= MAX_TUNNELS){
        is_finished = true;
      }
  //----
  int i = 0;    
  for(std::vector<Path*>::iterator iterator = paths.begin(); iterator != paths.end(); iterator++){
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
int path_optimization(Path* path, std::map<int, vertex*>& graph_points){
  //path = smoothing(path);
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
 
  smoothing(path);
  center_tunnel(path, graph_points);
  cut_tunnel(path, graph_points);

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  tunnel_optimization_time += duration_cast<microseconds>( t2 - t1 ).count();
  
  return center_tunnel_last_index_flag;
  //cut_tunnel(path);
  //path = center_tunnel(path);
}

void add_blocking_sphere(double* coordinates, double radius){
  static int blocking_spheres_counter = 0;
 vertex* blocking_sphere = new vertex(coordinates, 1, 0,radius, get_current_frame()); //don't care about index
 blocking_sphere->frame_index = get_current_frame();
 blocking_spheres_in_vertices.insert(std::pair<int,vertex*>(blocking_spheres_counter,blocking_sphere));
 blocking_spheres_counter++;
 rebuild_blocking_spheres_structure(coordinates, radius);
 

}

bool step_success_flag;
std::map<int, vertex*>& step(std::map<int, vertex*>& graph_points, vertex* new_vertex, vertex* parent, double goal_distance, 
  double collision_check_step, int parent_index, bool local_add_attempt){
 // print_vertex_coordinates(destination);
  static int counter = 0;
  
	int source_index = graph_points.size();
	double x_differential = new_vertex->location_coordinates[0] - parent->location_coordinates[0]; double kx = signum(x_differential); //program is using differential coordinates instead of absolute ones
	double y_differential = new_vertex->location_coordinates[1] - parent->location_coordinates[1]; double ky = signum(y_differential); 
	double z_differential = new_vertex->location_coordinates[2] - parent->location_coordinates[2]; double kz = signum(z_differential); 
    

	double l = compute_metric_eucleidean(new_vertex->location_coordinates,parent->location_coordinates, 3);        //distance of generated point(source) from nearest atom(destination)
	double l_planar = compute_metric_eucleidean(new_vertex->location_coordinates,parent->location_coordinates, 2); //distance of generated point from nearest atom in xy projection
	double l_from_source_to_goal = l - goal_distance;                                                               //distance by which the generated point is moved towards the nearest atom
  double alfa = atan(abs(z_differential) / l_planar);                                                             //angle of vertical and horizontal part of source to destiation trajectory
  double beta; x_differential != 0 ? beta = atan(abs(y_differential) / abs(x_differential)) : beta = 0;           //angle of planar parts 

  double new_x; double new_y; double new_z;


	if(!(l <= goal_distance)){
	 //geometric computation of new coordiantes which satisfy the maximal distance of generated point to nearest atom in rr tree
	 new_x = x_differential - kx * cos(alfa) * cos(beta) * l_from_source_to_goal;
	 new_y = y_differential - ky * cos(alfa) * sin(beta) * l_from_source_to_goal;
	 new_z = z_differential - kz * sin(alfa) * l_from_source_to_goal;
 
	 new_vertex->location_coordinates[0] = new_x + parent->location_coordinates[0]; new_vertex->location_coordinates[1] = new_y + parent->location_coordinates[1]; new_vertex->location_coordinates[2] = new_z + parent->location_coordinates[2];
  }
  if(is_in_obstacle(new_vertex->location_coordinates, probe_radius, CHECK_WITH_BLOCKING_SPHERES)){
    //static int deleted = 1;
    //cout << "del " << deleted << endl;
    //deleted++;
    delete new_vertex;
    step_success_flag = false;    
    return graph_points; //discard colliding point
  }
  
  parent->children_indices.push_back(new_vertex->index); //adds index of 
	graph_points.insert(std::pair<int,vertex*>(tree_index, new_vertex));
  
  step_success_flag = true;
  if(local_add_attempt && get_current_frame > 1){add_to_tree(new_vertex->location_coordinates, tree_index, local_priority_kdTree);} //local and global trees would be literally same in first
  add_to_tree(new_vertex->location_coordinates, tree_index, global_kdTree); //new node is always put in global tree
  tree_index++;

	if(termination_check(new_vertex->location_coordinates)){

    
		found = true;
    //new vertex is the last node added to terminated path
		Path* path;
    path = backtrack(graph_points, new_vertex->index);
    if(path == NULL){
      return graph_points;
    }
    
    int endpoint_node_index = path_optimization(path, graph_points);
   
    int is_duplicate = is_duplicated(path);

    bool add_sphere;
    //cout << "endp " << endpoint_node_index << endl;
    /*
    for(std::map<int, vertex*>::iterator it = path.begin(); it != path.end();it++){
      cout << it->second->index << endl;
    }
    */

    /*
    double* location_coordinates = new double[3];location_coordinates[0] = path->get_node_coordinates(path->ending_index)[0];location_coordinates[1] = path->get_node_coordinates(path->ending_index)[1]; location_coordinates[2] = path->get_node_coordinates(path->ending_index)[2];
    double blocking_sphere_radius = path->get_node_pointer(path->ending_index)->radius;
    */
    std::map<int, vertex*>& path_vertices = path->get_vertices();
    double* direction_vector = add_vectors(path_vertices[path->ending_index]->location_coordinates, path_vertices[path->get_parent_index(path->ending_index)]->location_coordinates, SUBTRACTION);
    double* centered_node_info = center_node(path->get_node_pointer(path->ending_index), NULL, NULL, direction_vector, test_sphere_radius, false, false, true); delete [] direction_vector;
    double* location_coordinates = new double[3]; location_coordinates[0] = centered_node_info[0]; location_coordinates[1] = centered_node_info[1]; location_coordinates[2] = centered_node_info[2];
    double blocking_sphere_radius = centered_node_info[3];
    delete [] centered_node_info;


    if(is_duplicate == -1){
      paths.push_back(path);
      duplicate_count.push_back(0);
      add_sphere = true;
    } else{
      if(!are_duplicated_pdbs_created){
        duplicate_count[is_duplicate] += 1;
        delete path;
        add_sphere = true;
      }
      else {
        duplicated_paths.push_back(path);
        duplicate_count[is_duplicate] += 1;
        add_sphere = true;
      }
    }
    //cout << "size " << location_coordinates.size();
    //cout << is_duplicate << endl;
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
    
    //graph_points = reset(graph_points);
    counter++;
    
	}
   
	return graph_points;
}

void print_vertex_coordinates(vertex* vert){
	std::cout <<"printing vertex location coordinates: x: " << vert->location_coordinates[0] << " y: " << vert->location_coordinates[1] << " z : " << vert->location_coordinates[2] << endl 
	<< "parent_index: " << vert->parent_index << endl;
}


MPNN::MultiANN<double>* build_local_tree(){
  int *topology = new int[dimension];
    MPNN::ANNcoord *scale = new MPNN::ANNcoord[dimension];
    for(int i=0;i<dimension;i++) {
      scale[i] = 1;
      topology[i] = 1;
     }
    MPNN::MultiANN<double>* local_kdTree = new MPNN::MultiANN<double>(dimension,1,topology,(MPNN::ANNpoint)scale);
    delete [] topology;
    delete [] scale;
    return local_kdTree;
};

MPNN::MultiANN<double>* new_frame_initalisation(map<int, vertex*>& graph_points){
  run_next_frame(); //switches frame number and loads next coordinate structure
  //checks which vertices from the old frame are not colliding in current frame, marks them as valid ones
  for (map<int, vertex*>::iterator iterator = graph_points.begin(); iterator != graph_points.end(); iterator++){
    if(!is_in_obstacle(iterator->second->location_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)){
      iterator->second->add_valid_frame(get_current_frame());
    }
  }
  int i = 1;
  //checks endpoints of paths, if they are valid in current frame. If succeeds, the path is automatically valid in this frame too, cause it
  //always possible to connect it to previous' frame's one
  for (vector<Path*>::iterator iterator = paths.begin(); iterator != paths.end(); iterator++){
    if(!is_in_obstacle((*iterator)->get_node_coordinates((*iterator)->ending_index), test_sphere_radius, DONT_CHECK_WITH_BLOCKING_SPHERES)){
      (*iterator)->add_valid_frame(get_current_frame());
    } else {
      cout << "path " << i << " is invalid in frame " << get_current_frame() << endl;
    }
    i++;
  }
  //first frame doesn't have priority tree, therefore it is not needed to delete it in second frame
  if(get_current_frame() > 2){
    delete local_priority_kdTree;
  }
  MPNN::MultiANN<double>* cur_frame_priority_kdTree = build_local_tree();

}


void build_tree(){
  
    int *topology = new int[dimension];
    MPNN::ANNcoord *scale = new MPNN::ANNcoord[dimension];
    for(int i=0;i<dimension;i++) {
      scale[i] = 1;
      topology[i] = 1;
 
    }
    global_kdTree = new MPNN::MultiANN<double>(dimension,1,topology,(MPNN::ANNpoint)scale);
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



void write_to_pdb(ofstream& file, std::map<int, vertex*>& atoms){
  std::map<int, vertex*>::iterator iterator;
  for(iterator = atoms.begin(); iterator != atoms.end(); iterator++){
      file << setprecision(3) << fixed;
      file << "ATOM";

      file.width(7); file << iterator->second->frame_index;
      float c1 = 1;
      file << "  N   ILE E  16";
      file.width(12); file <<  iterator->second->location_coordinates[0]; file.width(8); file <<  iterator->second->location_coordinates[1]; file.width(8); file <<  iterator->second->location_coordinates[2]; file.width(2); file.precision(2); file << "  "; file.width(4); file << c1;
      file.width(1); file << " "; file.width(5); file << iterator->second->radius; file.width(12); file << "N"<<endl;      
  }
}


void write_paths_to_pdbs(std::string filename, std::vector<Path*>& paths){

  int counter = 0;
  std::vector<Path*>::iterator iterator;
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

  std::map <int,vertex*> graph_points; //points cannot change the way they are accessed when others are added or deleted from main data structure, meaning vector representation would cause trouble, mapping is easiest way to index them
  load_parameters("config.txt");
  init_protein_struct();
  init_random();
  if(!is_init_position_given()){
    find_initial_position(0.3);
  }
  is_init_position_given() ? cout  << "init position from config file: " << qix << " " << qiy << " " << qiz << endl:cout  << "valid init position found: " << qix << " " << qiy << " " << qiz << endl;
  cout << get_current_frame() << endl;
  cout << "probe radius: " << probe_radius << endl;
  graph_points = init_start(graph_points);
  

  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  for(int i = 0; i < get_frames_count(); i++){
    if(i > 0) {
      new_frame_initalisation(graph_points);
    }
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    print_frame_borders();
  //---main part
    graph_points = find_path(graph_points, &num_iterations_stats_normal);
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

  std::ofstream blocking_file;
  blocking_file.open("blocking_spheres.pdb");
  write_to_pdb(blocking_file, blocking_spheres_in_vertices);
  blocking_file.close();
  
  
  int tree_size_at_finish = graph_points.size();
  clear_graph_points(graph_points);
  delete kdTree;
  graph_points.clear();

  
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

  string data = sout.str();

  cout << data;

  std::ofstream outfile;

  outfile.open("runtime_stats.txt", std::ios_base::app);
  outfile << data;
  outfile << endl;
  
 
  return 0;
}

