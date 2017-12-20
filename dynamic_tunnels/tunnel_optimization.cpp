#include <chrono>
#include <time.h>

#include "world.h"
#include "obstacle.h"
#include "support_methods.h"
#include "rrt_methods.h"
#include "MPNN/include/DNN/multiann.h"
#include "MPNN/include/DNN/ANN.h"

using namespace std;
using namespace std::chrono;




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

  cout << "path endpoint index " << path->get_endpoint_index() <<endl;
  cout << "deleted node index " << vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1] <<endl;
  int deleted_index = -1;
  if(path->get_endpoint_index() != vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1]){ //endpoint doesn't have child, would crash program
    deleted_index = path->get_child_index(vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1]);
  } else {
    return;
  }
  shared_ptr<vertex> cur = path->get_node_pointer(vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1]);

  int prev = path->get_size();
  cut_subtree_in_main_tree(path->get_node_pointer(deleted_index),false); //from element before the last
  path->erase_from_index(deleted_index);
  //cout << "path post " << prev - path->get_size() << endl;
  //exit(0);
  
}

void smooth_tunnel(shared_ptr<Path> path){
  path_nodes = path->get_vertices();
  shared_ptr<vertex> cur = path_nodes
  while(true){
    while(true){

    }
  } 
  
}

Function OptimizationPhaseOneAlgorithmOne(initialTunnel)
foreach i ∈ { 0, 1, ..., initialPath.size () − 1 } do
foreach j ∈ { 0, 1, ..., initialPath.size () − 1 } do
if collisionFree ( initialPath [ i ] , initialPath [ j ]) then
initialPath ← deleteNodes ( initialPath, i, j ) ;
end
end
end
return initialPath;

void center_tunnel(shared_ptr<Path> path){
  std::map<int, shared_ptr<vertex>>::iterator iterator;

  int counter = 0;
  int index = -1;
 
  std::map<int, shared_ptr<vertex>>& path_vertices = path->get_vertices();
  //cout << "center_tunnel: path_vertices size " << path_vertices.size() << endl;
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

      if(cur_node->get_children_count() > 0){
        index = cur_node->get_child_index();
      }
      continue;
    }


    shared_ptr<vertex> prev_node = path_vertices[cur_node->get_parent_index()];
    //cout << "center tunnel: counter " << counter << " cur_node " << cur_node->get_index() << endl;
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

