#include <chrono>
#include <time.h>
#include <queue>

#include "world.h"
#include "obstacle.h"
#include "support_methods.h"
#include "rrt_methods.h"
#include "MPNN/include/DNN/multiann.h"
#include "MPNN/include/DNN/ANN.h"

using namespace std;
using namespace std::chrono;

double SMOOTHING_SUCCESS_THRESHOLD = 5;
double MAXIMUM_SMOOTHING_DISTANCE = 0.8;






/*
Function which returns node approximately centered around the tunnel path. It can be used in two modes, either it centers node without respect to neighboring nodes, which leads to jagged path
but is more accurate or it takes into account two neighboring nodes and moving them in the same direction as the main node aimed for centering, which results in smoother path but more approximate
centering. Function returns new node without deleting the original, so the user must delete the original one if that is desired. Neighboring nodes' coordinates are changed only.
Function firslty finds plane perpendicular to supplied direction vector, using two perpendicular vectors to it and then performs "gradient descent" trying to move the node in the plane in attempt
to find the maximum possible radius to fit inside the tunnel. 
*/
double* center_node(shared_ptr<vertex> centered_node, shared_ptr<vertex> previous_node, shared_ptr<vertex> next_node, double* coordinate_vector, double increment, bool are_neighbors_centered, bool rewrite_node, bool return_node_info){
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  //sanity check
  if(isnan(centered_node->get_location_coordinates()[0])){
    cout << "NAN BEFORE CENTERING" << endl;
    create_segfault();
  }
  


  double normalised_length = 1;
  //what must be deleted perpendicular vector, node_coordinates, cur_perpendicular_vector, cur_second_perpendicular_vector, shift_vector, shifted_coordinates
  //first perpendicular vector, found by transposition of coordinates
 // cout << "coordinate vector " << endl;
  //print_vector(coordinate_vector);
  double perpendicular_vector_2D1[3]  = {-coordinate_vector[1], coordinate_vector[0], 0};
  double perpendicular_vector_2D2[3]  = {0, coordinate_vector[2], -coordinate_vector[1]};
  double* perpendicular_vector = add_vectors(perpendicular_vector_2D1, perpendicular_vector_2D2, ADDITION);
  //cout << "perpendicular_vector_2D1 " << endl;
  //print_vector(perpendicular_vector_2D1);
  //cout << "perpendicular_vector_2D2 " << endl;
  //print_vector(perpendicular_vector_2D2);
  //cout << "perpendicular_vector " << endl;
  //print_vector(perpendicular_vector);
  
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
    if(isnan(previous_node->get_location_coordinates()[0])){
        cout << "NAN BEFORE CENTERING, PREV" << endl;
        create_segfault();
      }
    if(next_node != NULL){
      next_node_shifted_coordinates = next_node->copy_coordinates(); //happens for last node, that one should be centered, unlike the root node
      if(isnan(next_node->get_location_coordinates()[0])){
        cout << "NAN BEFORE CENTERING, NEXT" << endl;
        create_segfault();
      }
    }
      
  } 
  
  //maximum radius in default position, used in first iteration to find larger one 
  //cout << "first " << endl;

  
  double cur_radius = approximate_radius_custom_frame(node_coordinates, centered_node->get_radius(), 0.01, 0.1, centered_node->get_first_frame());
  //sanity check
  if(get_current_frame() == 2){
    //cout << "index " << centered_node->get_index() << endl;
    //cout << centered_node->get_frame_index() << endl;
    //cout << "first " << approximate_radius_custom_frame(node_coordinates, centered_node->get_radius(), 0.01, 0.1, 1) << endl;
    //cout << "second " << approximate_radius_custom_frame(node_coordinates, centered_node->get_radius(), 0.01, 0.1, 2) << endl;
  //sanity check
  //sanity check
  }

  if(cur_radius == 0){
    cout << "centering node which is already colliding!" << endl;
    create_segfault();
  }
  //sanity check
  if(is_in_obstacle_custom_frame(node_coordinates, cur_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, centered_node->get_first_frame())){
    cout << "cur_radius_colliding pre" << endl;
    create_segfault();
  }
  //for comparison
  double check_radius = cur_radius;
  double last_valid_radius = cur_radius;

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

      //... and maximum radius of ball in that placement is found, nodes with smaller than the largest last radius should be automatically discarded here
      cur_radius = approximate_radius_custom_frame(shifted_coordinates, last_valid_radius, 0.01, 0.1, centered_node->get_first_frame());
      //node is colliding
      if(cur_radius == 0){
        delete [] cur_perpendicular_vector; delete [] cur_second_perpendicular_vector; delete [] shift_vector; delete [] shifted_coordinates; 
        continue;
      }
      
      //sanity check
      if(is_in_obstacle_custom_frame(shifted_coordinates, cur_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, centered_node->get_first_frame())){
        cout << "cur_radius_colliding" <<endl;
        create_segfault();
      }
      //maximum radius in this position is checked subtracted from maximum radius found in previous iteration and checked against the maximum increase found in current iteration, which starts at zero,
      //therefore the new maximum radius must be greater than the previous one
      if((cur_radius - check_radius) > max_increase && !is_in_obstacle_custom_frame(shifted_coordinates, cur_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, centered_node->get_first_frame())) {
        
        //cout << "sanity check " << are_neighbors_centered << endl;
        //--- can be buggy, check later
        if(are_neighbors_centered){
          //new coordinates shifted in the same way as main node, this is to prevent excessive deviation and resulting jaggedness of tunnel
          double* neighbors_shift_vector = copy_vector(shift_vector); multiply_vector(neighbors_shift_vector, 0.5);
         // cout << neighbors_shift_vector[0] << endl;
         // cout << shift_vector[0] << endl;
          double* cur_prev_node_shifted_coordinates = add_vectors(neighbors_shift_vector, prev_node_shifted_coordinates, ADDITION);
          double* cur_next_node_shifted_coordinates = add_vectors(neighbors_shift_vector, next_node_shifted_coordinates, ADDITION);
          delete [] neighbors_shift_vector;
          
          //cout << "prev "<< is_in_obstacle(prev_node_shifted_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES);
          //print_vector(prev_node_shifted_coordinates);
          //cout << "next "<< is_in_obstacle(next_node_shifted_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES);
          //print_vector(next_node_shifted_coordinates);
          if((next_node == NULL && !is_in_obstacle_custom_frame(prev_node_shifted_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, previous_node->get_last_valid_frame()))||
            (next_node != NULL && !is_in_obstacle_custom_frame(prev_node_shifted_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, previous_node->get_last_valid_frame()) 
              && !is_in_obstacle_custom_frame(next_node_shifted_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, next_node->get_last_valid_frame()))){
             //test was succesfull, node can be moved
             //new radius is assigned
             max_increase = cur_radius - check_radius;
             last_valid_radius = cur_radius;
             delete [] node_coordinates;
             node_coordinates = copy_vector(shifted_coordinates);
             //cout << "node coordinates: " << endl;
             //print_vector(node_coordinates);
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
          max_increase = cur_radius - check_radius;
          last_valid_radius = cur_radius;
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

  //sanity check
  if(last_valid_radius < check_radius){
    cout << "entering loop assigned smaller radius than in beginning!" << endl;
    create_segfault();
  }
  double ret_radius;
  last_valid_radius > check_radius ? ret_radius = last_valid_radius : ret_radius = check_radius;
  delete [] perpendicular_vector;
  if(rewrite_node){
   // cout << "rewriting " << endl;
   // cout << "radius " << last_valid_radius << endl;
   // print_vector(node_coordinates);
  centered_node->set_coordinates(node_coordinates);
  centered_node->set_radius(ret_radius);
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
    if(is_in_obstacle_custom_frame(centered_node->get_location_coordinates(), centered_node->get_radius(), DONT_CHECK_WITH_BLOCKING_SPHERES, centered_node->get_first_frame())){
      cout << "centered node colliding!" <<endl;
      create_segfault();
    }
    if(isnan(centered_node->get_location_coordinates()[0])){
      cout << "NAN after CENTERING" << endl;
      create_segfault();
    }
    return NULL;
  } else {
    double* ret = new double[4];
    ret[0] = node_coordinates[0]; ret[1] = node_coordinates[1]; ret[2] = node_coordinates[2]; ret[3] = cur_radius;
    delete [] node_coordinates;
    if(is_in_obstacle_custom_frame(centered_node->get_location_coordinates(), centered_node->get_radius(), DONT_CHECK_WITH_BLOCKING_SPHERES, centered_node->get_first_frame())){
      cout << "centered node colliding!" <<endl;
      create_segfault();
    } 
    if(isnan(centered_node->get_location_coordinates()[0])){
      cout << "NAN BEFORE CENTERING" << endl;
      create_segfault();
    }
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

int is_node_duplicate(shared_ptr<vertex> cur_node, bool is_local){

  int exists_in_path; 
  if(is_local){
    exists_in_path = local_priority_kdTree_coordinates->get_node_pointer(cur_node->get_index())->is_in_path();
  } else {
    exists_in_path = global_tree_points->get_node_pointer(cur_node->get_index())->is_in_path();
  }

  if(exists_in_path != -1){ 
    if (cur_node->get_first_frame() == paths[exists_in_path]->get_node_pointer(cur_node->get_index())->get_first_frame()){
      return exists_in_path;
    } else {
      return -1;
    }
  } else {
    return -1;
  }

}


int cut_tunnel(shared_ptr<Path> path){
  std::map<int, shared_ptr<vertex>>& path_vertices = path->get_vertices();
  std::vector<int> vertex_indices_to_be_deleted;

  int endpoint_index = path->get_endpoint_index();
  int endpoint_parent_index = path->get_parent_index(endpoint_index);
  bool is_local = path_vertices[endpoint_index]->is_local();
  //print_map_indices(path->get_vertices(), "cut tunnel");
  while((endpoint_index != -1 && endpoint_parent_index != -1) && is_node_duplicate(path_vertices[endpoint_index], is_local) == -1){   
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
    double* node_information = center_node(path_vertices[endpoint_index], NULL, NULL, direction_vector, test_sphere_radius / 2, false, false, true);
    if(node_information [3] > test_sphere_radius){      
      vertex_indices_to_be_deleted.push_back(endpoint_index);
    } 

    delete [] direction_vector;
    delete [] node_information;

    endpoint_index = path->get_parent_index(endpoint_index);
    endpoint_parent_index = path->get_parent_index(endpoint_index);

  }

  //cout << "path endpoint index " << path->get_endpoint_index() <<endl;
  //cout << "deleted node index " << vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1] <<endl;
  int deleted_index = -1;
  if(vertex_indices_to_be_deleted.size() == 0){
    return -2;
  }
  if(path->get_endpoint_index() != vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1]){ //endpoint doesn't have child, would crash program
    deleted_index = path->get_child_index(vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1]);
  } else {
    return -2;
  }
  shared_ptr<vertex> cur = path->get_node_pointer(vertex_indices_to_be_deleted[vertex_indices_to_be_deleted.size() - 1]);

  int prev = path->get_size();
  //cut_subtree_in_main_tree(path->get_node_pointer(deleted_index),false); //from element before the last
  path->erase_from_index(deleted_index);
  return deleted_index;
  //cout << "path post " << prev - path->get_size() << endl;
  //exit(0);
  
}

bool search(shared_ptr<vertex> base, shared_ptr<vertex> target){

    int frame_diference = target->get_first_frame() - base->get_first_frame();
    double  direction_vector_length = compute_metric_eucleidean(target->get_location_coordinates(), base->get_location_coordinates());
    int num_of_interpolated_steps = direction_vector_length / min_step;

    if(direction_vector_length <= min_step){
        return false;
    }

    double* increment_vector = create_direction_vector(base, target, min_step);

    vector<vector<int>> matrix;
    for(int i = 0; i < frame_diference; i++){
      vector<int> matrix_row;
      matrix.push_back(matrix_row);
      for(int j = 0; j < num_of_interpolated_steps; j++){
        matrix[i].push_back(-1);
      }
    }
    cout << "gone"  << endl;
    for(int i = 0; i < num_of_interpolated_steps; i++){
        double* current_coordinates = add_vectors(base->get_location_coordinates(), increment_vector, ADDITION);
        for(int j = 0; j < frame_diference; j++){
            if(is_in_obstacle_custom_frame(current_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, base->get_first_frame() + j)){
                matrix[j][i] = 0;
            }
            else {
                matrix[j][i] = 1;
            } 
        }
        delete [] current_coordinates;
    }
   
    queue <array<int, 2>> visited;
    array<int, 2> cur = {0, 0};
    visited.push(cur);

    bool success = false;

    int counter = 0;

    while(visited.size() != 0){
        if(cur[0] < frame_diference - 1){
            if(matrix[cur[0] + 1][cur[1]] == 1){
                if(cur[0] == num_of_interpolated_steps - 1 && cur[1] == frame_diference - 1){
                    success = true;
                    break;
                }
                array<int, 2> new_position = {cur[0] + 1, cur[1]};
                visited.push(new_position);
            }
        }   

        if(cur[0] < num_of_interpolated_steps - 1){
            if(matrix[cur[0]][cur[1] + 1] == 1){
                if(cur[0] == num_of_interpolated_steps - 1 && cur[1] == frame_diference - 1){
                    success = true;
                    break;
                }
                array<int, 2> new_position = {cur[0], cur[1] + 1};
                visited.push(new_position);
            }
        }

    cur = visited.front();
    visited.pop();
    counter++;
    if(counter > 10000000){
      cout << "STUCK" << endl;
      create_segfault();
    }

    }

    cout << "gone2"  << endl;
    delete [] increment_vector;
    return success;
}



bool delete_nodes_inbetween(shared_ptr<Path> path, shared_ptr<vertex> base, shared_ptr<vertex> target){
  //means that there was time travel inbetwee

  double  direction_vector_length = compute_metric_eucleidean(target->get_location_coordinates(), base->get_location_coordinates());
 // cout << "direction vector length " << direction_vector_length << endl;
  if(direction_vector_length <= min_step){
    return false;
  } 
  double* increment_vector = create_direction_vector(base, target, min_step);

  //shared_ptr<vertex> cur = first;
 // vector<shared_ptr<vertex>> interpolated_nodes;

  //check if straight line is colliding
  //applying to new found path nodes, no need to rebuild kdTrees and all new found path nodes are local
  //until the length to target is smaller than min_step, try to add new node to straight line trajectory between and check for collisions
  double* new_position = base->get_location_coordinates();
  double* new_position_prev = copy_vector(new_position);
  while(compute_metric_eucleidean(target->get_location_coordinates(), new_position_prev) > min_step){
    new_position = add_vectors(new_position_prev, increment_vector, ADDITION);
    delete [] new_position_prev;
    if(!is_in_obstacle_custom_frame(new_position, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, base->get_first_frame())){
 //     cur = make_shared<vertex>(new_position, 0, probe_radius, get_current_frame(), true);
 //     interpolated_nodes.push_back(cur);
 //     cout << "distance " << compute_metric_eucleidean(target->get_location_coordinates(), new_position) <<endl;
      new_position_prev = copy_vector(new_position);
      delete [] new_position;
    } else {
//      cout << "colliding" << endl;
      delete [] new_position;
      delete [] increment_vector;
      return false;
    }
  }

  //..if not, delete all nodes in between 
  shared_ptr<vertex> cur = base->get_child_pointer().lock();

//  cout << "first " << cur->get_index() << endl;
//  cout << "second " << target->get_index() << endl;
//  cout << "endpoint " << path->get_endpoint_index() << endl;


  while(cur->get_index() != target->get_index()){
    path->erase_node_with_reconnecting(cur->get_index());
//    cout << "in cur " << cur->get_index() << endl;
 //   cout << "in second " << target->get_index() << endl;
 //   cout << "in child " << cur->get_child_index() << endl;
    cur = cur->get_child_pointer().lock();  
  }
  
  delete [] increment_vector;
  delete [] new_position_prev;
  return true;
}

bool delete_nodes_inbetween_single(shared_ptr<Path> path, shared_ptr<vertex> base, shared_ptr<vertex> middle, shared_ptr<vertex> target){
  if(base->get_first_frame() < target->get_first_frame()){
    //cout << "IM HERE" << endl;
    if(search(base, target)){
      path->erase_node_with_reconnecting(middle->get_index());
    }
  } else {
      delete_nodes_inbetween(path, base, target);
  } 
}

double smoothing_vol2(shared_ptr<Path> path){
  int beginning_index = path->get_beginning_index(); 
  //cout << path.size() << endl;
  if(path->get_size() < 3){
    return 0;
  }

  shared_ptr<vertex> cur = path->get_beginning_node();
  shared_ptr<vertex> child = cur->get_child_pointer().lock()->get_child_pointer().lock();
  int success_counter = 0;
  int run_counter = 0;
  while(true){

    //too much space inbetween can lead trajectory which is not described by the spheres, there would be just empty space inbetween, which doesn't really matter
    //for the trajectory itself, becaouse important information are extracted from it before the procedure, but is bad for visualization. It's up user to edit this parameter.
   //cout <<  "dist " << compute_metric_eucleidean(cur->get_location_coordinates(), child->get_location_coordinates()) <<endl;
   // cout <<  "param " << MAXIMUM_SMOOTHING_DISTANCE << endl;
    if(compute_metric_eucleidean(cur->get_location_coordinates(), child->get_location_coordinates()) < MAXIMUM_SMOOTHING_DISTANCE){
      if(delete_nodes_inbetween_single(path, cur, cur->get_child_pointer().lock(), child)){
        success_counter++; 
      }
    }  
    run_counter++; 
    
    cur = child;
    if(child->get_child_pointer_null_permisive().lock() == NULL || child->get_child_pointer().lock()->get_child_pointer_null_permisive().lock() == NULL){
      break;
    }
    child = child->get_child_pointer().lock()->get_child_pointer().lock();
  }
  
  //cout << ((double) success_counter / (double) run_counter) * 100 << endl;
  return  ((double) success_counter / (double) run_counter) * 100;
}


void add_interpolated_node(shared_ptr<Path> path, shared_ptr<vertex> parent, shared_ptr<vertex> cur, shared_ptr<vertex> child){
    int index = tree_index;
    tree_index++;
    cur->set_index(index);
    cur->set_child_pointer(child);
    cur->set_parent_pointer(parent);
    path->add_node(cur);
}

/*
void interpolate_nodes(shared_ptr<Path> path, shared_ptr<vertex> base, shared_ptr<vertex> target){

  double  direction_vector_length = compute_metric_eucleidean(base->get_location_coordinates(), target->get_location_coordinates());
  //there isn't any place to for a new node, happens if there wasn't any reconnection
  if(direction_vector_length <= min_step){
    return;
  }
 // cout << "OLD PATH SIZE " << path->get_vertices().size() << endl; 
  //the increment vector among which the path segment will be interpolated
  double* increment_vector = add_vectors(target->get_location_coordinates(), base->get_location_coordinates(), SUBTRACTION);  
  normalize_vector(increment_vector, min_step);

  //check special case when there's only one new node to add
  bool is_singleton;
  direction_vector_length <= 2*min_step ? is_singleton = true : is_singleton = false;
  shared_ptr<vertex> parent = base;
  
  //adding new position for cur node
  double* new_position = add_vectors(base->get_location_coordinates(), increment_vector, ADDITION);
  shared_ptr<vertex> cur =  make_shared<vertex>(new_position, 0, probe_radius, get_current_frame(), true);

  //case when no node is added is already resolved and adding the first child is nonstandard operation happening in every other relevant case
  base->set_child_pointer(cur);
  
  //singleton case
  if(is_singleton){
    target->set_parent_pointer(cur);
    add_interpolated_node(path, base, cur, target);
    delete [] new_position;
    return;
  }

  //initiating the cur's child vertex, so pointer to it can be added in current loop cycle
  double* new_position_child = add_vectors(new_position, increment_vector, ADDITION);
  shared_ptr<vertex> child =  make_shared<vertex>(new_position_child, 0, probe_radius, get_current_frame(), true);
  delete [] new_position_child;
  delete [] new_position;


  //while the child's distance to target node is longer than min_step
  int counter = 0;
 // print_vector(increment_vector);
  while(compute_metric_eucleidean(child->get_location_coordinates(), target->get_location_coordinates()) > min_step){

   //   print_vector(child->get_location_coordinates());
   //   print_vector(target->get_location_coordinates());
      add_interpolated_node(path, parent, cur, child);

      parent = cur;
      cur = child;
      //position of new child node
      new_position_child = add_vectors(cur->get_location_coordinates(), increment_vector, ADDITION);
      child = make_shared<vertex>(new_position_child, 0, probe_radius, get_current_frame(), true);
      //position is now copied in child vertex
      delete [] new_position_child;
  }
  //final connection
  add_interpolated_node(path, parent, cur, child);
  add_interpolated_node(path, cur, child, target);
  target->set_parent_pointer(child);
  if(!test_path_noncolliding_static(path)){
//    cout << "INTERPOLATED PATH COLLIDING " << endl;
//    exit(0);
  }
//  cout << "NEW PATH SIZE " << path->get_vertices().size() << endl;
//  cout << "STUFF " << child->get_index() << " " << target->get_index() << endl; 
}
*/
/*
void interpolate_smoothed_segments(shared_ptr<Path> path){
  if(path->get_size() < 2){
//    cout << "TEST ERROR: interpolate_smoothed_segments: path has less than two nodes! " << endl;
    return;
  }

  shared_ptr<vertex> base = path->get_beginning_node();
  shared_ptr<vertex> target;

//  cout << "ENDPOINT " << path->get_endpoint_index() << endl;
  while(base->get_index() != path->get_endpoint_index()){
    target = base->get_child_pointer().lock();
    interpolate_nodes(path, base, target);
    base = target;

//    cout << "BASE " << base->get_index() << endl;
//    cout << "ENDPOINT NOW " << path->get_endpoint_index() << endl;
  }
}
*/


/*
void smooth_tunnel(shared_ptr<Path> path){
  shared_ptr<vertex> cur = path->get_beginning_node();
  
  int counter = 0;
  while(cur->get_index() != path->get_endpoint_index()){
    if(cur->is_in_some_valid_path()){
      cur = cur->get_child_pointer().lock();
      continue;
    }
    shared_ptr<vertex> cur_inner = path->get_beginning_node();
    bool is_inner_first = true;
    while(cur_inner != NULL){
//      cout << endl << "cur " << cur->get_index() << endl;
//      cout << "cur_inner " << cur_inner->get_index() << endl;
//      cout << "distance " << compute_metric_eucleidean(cur->get_location_coordinates(), cur_inner->get_location_coordinates()) << endl;

      if(is_inner_first){
        if(cur->get_index() == cur_inner->get_index()){
          is_inner_first = false;
          cur_inner = cur_inner->get_child_pointer().lock();
          continue;
        } else {
          delete_nodes_inbetween(path, cur_inner, cur);
        }
      }
      else {
      delete_nodes_inbetween(path, cur, cur_inner);
      }
      cur_inner = cur_inner->get_child_pointer_null_permisive().lock();
    }
    
    cur = cur->get_child_pointer().lock();
  }

  interpolate_smoothed_segments(path);
  
}
*/




int center_tunnel(shared_ptr<Path> path){
  //std::map<int, shared_ptr<vertex>>::iterator iterator;

  int counter = 0;
  int index = -1;
 
  //std::map<int, shared_ptr<vertex>>& path_vertices = path->get_vertices();
  // cout << "size " << path_vertices.size() << endl;
  //cout << "center_tunnel: path_vertices size " << path_vertices.size() << endl;
  if(path->get_size() < 2) return -2;



  shared_ptr<vertex> cur_node = path->get_beginning_node()->get_child_pointer().lock();
  double* direction_vector = add_vectors(path->get_beginning_node()->get_location_coordinates(),cur_node->get_location_coordinates(), SUBTRACTION);
  center_node(path->get_beginning_node(), NULL, NULL, direction_vector, 0.5, false, true, false);
  delete [] direction_vector;


  //shared_ptr<vertex> cur_node = path->get_endpoint_node();

  //path->print_path();
  while(true){
    if(get_current_frame() == 2){
      //cout <<"shit "<< cur_node->get_index() << endl; 
    }
    
    /*
    With reused tree, that already centered node is inside some previous trajectory. Therefore instead of costy centering procedure it can just be copied from that path.
    exists_in_path() method returns -1 if the node is new. else it returns index into path in paths vector.
    Important - node must be valide in the same frame as current node, because molecular positions are different in every frame!
    Duplicated paths are deleted, therefore we must be sure to not index into them, this is achieved by calling the "labeling" function, which points the nodes in main data structures
    into valid paths, when the path is added into main paths vector 
    */
    //check if node exists in some other path and in the same frame
    int related_node_path_index = is_node_duplicate(cur_node, cur_node->is_local());

    if(related_node_path_index != -1){
      if(get_current_frame() == 2){
        //cout <<"shit2 "<< cur_node->get_index() << endl; 
      }
      shared_ptr<vertex> related_node = paths[related_node_path_index]->get_node_pointer(cur_node->get_index());
      
      cur_node->set_coordinates(related_node->get_location_coordinates()); //copy identical nodes coordinates to current;
      cur_node->set_radius(related_node->get_radius());

      //cout << cur_node->get_index() << endl;

      if(cur_node->get_children_count() > 0){
        if(path->is_node_endpoint(cur_node)){
          return -2;
        }
        cur_node = cur_node->get_child_pointer_null_permisive().lock();
        if(cur_node == NULL){
          cout << "TEST ERROR, center tunnel returns empty node:571 " <<endl;
          create_segfault();
        }
      }
      continue;
    }


    shared_ptr<vertex> parent_node = cur_node->get_parent_pointer().lock();
    shared_ptr<vertex> child_node = NULL;
    if(!path->is_node_endpoint(cur_node)){
      child_node = cur_node->get_child_pointer_null_permisive().lock();
      if(child_node == NULL){
        cout << "TEST ERROR, center tunnel returns empty node " <<endl;
        create_segfault();
      }
    }

    double* direction_vector = add_vectors(cur_node->get_location_coordinates(), parent_node->get_location_coordinates(), SUBTRACTION);
    if(vector_length(direction_vector) == 0){
      create_segfault();
    }
    

    
    if(child_node != NULL){
      
      if(isnan(cur_node->get_location_coordinates()[0])){
        create_segfault();
      }
      center_node(cur_node, parent_node, child_node, direction_vector, 0.5, false, true, false);
      if(isnan(cur_node->get_location_coordinates()[0])){
        create_segfault();
      }


    } else {
      if(isnan(cur_node->get_location_coordinates()[0])){
        create_segfault();
      }
      center_node(cur_node, parent_node , NULL, direction_vector, 0.5, false, true, false);
      if(isnan(cur_node->get_location_coordinates()[0])){
        create_segfault();
      }
    }

    delete [] direction_vector;
    
    /* deleting subtrees is important because of path blocking. Blocking sphere needs to be placed in the center of the tunnel at it's endpoint
    and the tree after this endpoint would create would create valid growth into the empty region, which obviously sucks. It's
    caused by the fact that the part of the tunnel growing after original endpoint ends after centering in empty space.
    both local and global coordinates of the tree after endpoint needs to be declared invalid now and deleted from kd trees,
    and the kd trees are rebuilded. In first frame, it's possible to only care about local structures, which speeds up static detection.
    */

    if(path->is_node_endpoint(cur_node)){
      return -2;
    }

    if(cur_node->get_radius() > test_sphere_radius){
      int child_index = cur_node->get_child_pointer().lock()->get_index();    
      if(child_index == -1){
          cerr << "your start position already lies in empty space!" << endl;
          return -2;      
      }
      //cut_subtree_in_main_tree(path->get_node_pointer(child_index), false); 
      path->erase_from_index(child_index);
      return child_index;
  }; //to_do

    
    if(cur_node->get_children_count() > 0){
      cur_node = cur_node->get_child_pointer().lock();
    } else {
      cout << "TEST ERROR, SOME SHIT HAPPEND " <<endl;
      create_segfault();
    }

  }
  
   
}

void center_tunnel_without_erasing(shared_ptr<Path> path){

  if(path->get_size() < 2) return;
  shared_ptr<vertex> cur = path->get_beginning_node()->get_child_pointer().lock();
  shared_ptr<vertex> parent = NULL;
  shared_ptr<vertex> child = NULL;

  //cout << "size " << path_vertices.size() << endl;
  //cout << "center_tunnel: path_vertices size " << path_vertices.size() << endl;
  while(true){  
    /*
    With reused tree, that already centered node is inside some previous trajectory. Therefore instead of costy centering procedure it can just be copied from that path.
    exists_in_path() method returns -1 if the node is new. else it returns index into path in paths vector.
    Important - node must be valide in the same frame as current node, because molecular positions are different in every frame!
    Duplicated paths are deleted, therefore we must be sure to not index into them, this is achieved by calling the "labeling" function, which points the nodes in main data structures
    into valid paths, when the path is added into main paths vector 
    */

    bool is_endpoint = false;
    parent = cur->get_parent_pointer().lock();
    if(cur->get_index() == path->get_endpoint_index()) is_endpoint = true;
    if(!is_endpoint) child = cur->get_child_pointer().lock();

    double* direction_vector = add_vectors(cur->get_location_coordinates(), parent->get_location_coordinates(), SUBTRACTION);
  
    //if(!is_endpoint) cout << "child_node " << child->get_index() << endl;
    center_node(cur, parent, child, direction_vector, 1, false, true, false);
    delete [] direction_vector;
    if(is_endpoint){
      break;
    }

    cur = child;
  }
}



int path_optimization(shared_ptr<Path> path){
  //path = smoothing(path);
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  //cout << "path size bef " << path->get_size() << endl;
  //smooth_tunnel(path);
  if(!test_path_noncolliding_static(path)){
      cout << "path colliding opt bef processing!" <<endl;
      create_segfault();
    }

  int delet_index_after_centering = center_tunnel(path);
  //cout << "path size aft_ctr " << path->get_size() << endl;
  if(!test_path_noncolliding_static(path)){
      cout << "path colliding opt_after_centr processing!" <<endl;
      create_segfault();
    }
  delet_index_after_centering = cut_tunnel(path);
  //cout << "path size aft_cut " << path->get_size() << endl;

  if(!test_path_noncolliding_static(path)){
      cout << "path colliding opt_post processing!" <<endl;
      create_segfault();
    }
  
  

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  tunnel_optimization_time += duration_cast<microseconds>( t2 - t1 ).count();
  
  return delet_index_after_centering;
  //cut_tunnel(path);
  //path = center_tunnel(path);
}


double find_tunnel_length_and_bottleneck(shared_ptr<Path> path){
  shared_ptr<vertex> prev = path->get_beginning_node();
  shared_ptr<vertex> cur = prev->get_child_pointer().lock();
  double length = 0;
  double bottleneck = DBL_MAX;

  while(cur != NULL){
    length += compute_metric_eucleidean(prev->get_location_coordinates(), cur->get_location_coordinates());
    if(cur->get_radius() < bottleneck){
      bottleneck = cur->get_radius();
    }
    prev = cur;
    cur = cur->get_child_pointer_null_permisive().lock(); 
  }
}



void path_optimization_postprocessing(shared_ptr<Path> path){
  //smooth_tunnel(path);
  //center_tunnel_without_erasing(path);


  set_path_bottleneck(path);
  if(!test_path_noncolliding_static(path)){
      cout << "path colliding pre post-processing!" <<endl;
      create_segfault();
  }
  double success = 100;
  do {
    //cout << "path_size pre" << path->get_size() << endl; 
    success = smoothing_vol2(path);
    //cout << "success " << success << endl;
    //cout << "path_size post" << path->get_size() << endl; 
  }
  while(success > SMOOTHING_SUCCESS_THRESHOLD);
  set_path_length(path);
  if(!test_path_noncolliding_static(path)){
      cout << "path colliding post post-processing!" <<endl;
      create_segfault();
  }

}

