#include <chrono>
#include <time.h>

#include "world.h"
#include "obstacle.h"
#include "support_methods.h"
#include "rrt_methods.h"


using namespace ozcollide;
using namespace std::chrono;
using namespace std;



 
//collision checking without checking the blocking spheres structure is used in functions not involving path search (path optimization and chcecking whether new node is inside the protein) 
const int CHECK_ONLY_BLOCKING_SPHERES = 2; 
const int CHECK_WITH_BLOCKING_SPHERES = 1; 
const int DONT_CHECK_WITH_BLOCKING_SPHERES = 0;


template <class T>
bool is_in_obstacle_basic(T loc_coord, double radius, int check_with_blocking_spheres, int frame){
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
      returned = get_frame_protein_tree(frame)->isCollideWithSphere(sphere); 
      t2 = high_resolution_clock::now(); 
      collision_check_time += duration_cast<microseconds>( t2 - t1 ).count(); 
      return returned; 
    break; 
    case CHECK_WITH_BLOCKING_SPHERES: 
     bool has_collided_with_protein_structure = get_frame_protein_tree(frame)->isCollideWithSphere(sphere); 
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

bool is_in_obstacle(double* loc_coord, double radius, int check_with_blocking_spheres){  
  return is_in_obstacle_basic(loc_coord, radius, check_with_blocking_spheres, get_current_frame());       
}

bool is_in_obstacle(std::array<double, 3> loc_coord, double radius, int check_with_blocking_spheres){ 
  return is_in_obstacle_basic(loc_coord, radius, check_with_blocking_spheres, get_current_frame());       
}

bool is_in_obstacle_custom_frame(double* loc_coord, double radius, int check_with_blocking_spheres, int frame){ 
  return is_in_obstacle_basic(loc_coord, radius, check_with_blocking_spheres, frame);      
}





double approximate_radius_basic(double* coordinates, double default_radius, double precision, double default_approximation_step, int frame){
  //cout << "frame " << frame << endl;

  double approximation_step = default_approximation_step;
  double cur_radius = default_radius + approximation_step;
  double last_valid_radius = default_radius;
  if(is_in_obstacle_custom_frame(coordinates, last_valid_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, frame)){
     return 0;
  }
  while(approximation_step > precision){
    double test_radius = cur_radius + approximation_step;
    bool hasCollided = is_in_obstacle_custom_frame(coordinates, test_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, frame);
    if(!hasCollided){
      cur_radius = test_radius;
      last_valid_radius = test_radius;
    }
    else{
      approximation_step = approximation_step / 2;
      cur_radius -= approximation_step;
      
    }
  }
  if(is_in_obstacle_custom_frame(coordinates, last_valid_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, frame)){
      cout << "approximated radius colliding!" <<endl;
      create_segfault();
    }
  return last_valid_radius;
}

//function which returns approximation of maximal radius of sphere fitting into the tunnel, expects non-colliding sphere 
double approximate_radius(double* coordinates, double default_radius, double precision, double default_approximation_step){
  return approximate_radius_basic(coordinates, default_radius, precision, default_approximation_step, get_current_frame());
}

//function which returns approximation of maximal radius of sphere fitting into the tunnel, expects non-colliding sphere 
double approximate_radius_custom_frame(double* coordinates, double default_radius, double precision, double default_approximation_step, int frame){
  return approximate_radius_basic(coordinates, default_radius, precision, default_approximation_step, frame);
}




int print_colliding(std::map<int ,shared_ptr<vertex>>& map){
  int counter = 0;
  for(std::map<int, shared_ptr<vertex>>::iterator iterator = map.begin(); iterator != map.end(); iterator++){
    if(is_in_obstacle(iterator->second->get_location_coordinates(), probe_radius, CHECK_ONLY_BLOCKING_SPHERES)){
    //    cout <<  " "  << iterator->second->get_index() << " ";
        counter++; 
    }
   }

   return counter;
}



void purge_nodes_in_blocking_spheres_recursion(shared_ptr<vertex> cur){
  int index = 0;
  //while there are unvisited children
  while(index < cur->get_children_count()){
    shared_ptr<vertex> next = cur->get_child_pointer(index).lock();
    bool is_next_purged = is_in_obstacle(next->get_location_coordinates(), probe_radius, CHECK_ONLY_BLOCKING_SPHERES);

    //next is in blocking sphere, therefore whole it's while branch will be deleted, otherwise we just move to next child
    if(is_next_purged) {
      //cout << "purged " << next-> get_index() << endl; 
      cut_subtree_in_main_tree(next, false, false);
    }
    else {
      purge_nodes_in_blocking_spheres_recursion(next); 
      index++;
    }
  }
}


void purge_nodes_in_blocking_spheres(shared_ptr<Path> optimizedTunnel){
  //first node after root in the tree structure
  shared_ptr<vertex> cur = local_priority_kdTree_coordinates->get_node_pointer(optimizedTunnel->get_beginning_node()->get_child_pointer().lock()->get_index());
  purge_nodes_in_blocking_spheres_recursion(cur);
  rebuild_kd_tree(true, true, true); 

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
    shared_ptr<vertex> next = cur->get_child_pointer(cur->get_children_count() -1).lock(); 
    cut_subtree_in_main_tree_recursion(next, has_debugging_output, next->is_local()); 
    if(has_debugging_output){ 
      cout << "aft "<< cur->get_children_count() << endl;  
    } 
    cur->delete_child(cur->get_children_count() -1); //last one remaining 
  } 
  if(has_debugging_output){cout << "deleting " <<cur->get_index()<<endl;} 
  local_priority_kdTree_coordinates->delete_node(cur->get_index()); 
 
} 
 
//deletes subtree beginning at index and deletes reference to subtree root from it's parent node  
void cut_subtree_in_main_tree(shared_ptr<vertex> first_delet_node, bool has_debugging_output, bool is_nearest_neighbor_search_structure_rebuilt){ 
  shared_ptr<vertex> node_in_tree; 
  shared_ptr<vertex> parent_node_in_tree; 
 
  node_in_tree = local_priority_kdTree_coordinates->get_node_pointer(first_delet_node->get_index()); 
  parent_node_in_tree = node_in_tree->get_parent_pointer().lock(); 
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
 
  int local_tree_size = local_priority_kdTree_coordinates->get_size(); 
  int prev = local_priority_kdTree_coordinates->get_size(); 
  cut_subtree_in_main_tree_recursion(node_in_tree, has_debugging_output, first_delet_node->is_local()); 

  if(is_nearest_neighbor_search_structure_rebuilt) rebuild_kd_tree(true, true, true); 
  //cout << "tree rebuilt" << endl; 
}

//Later I should research how to use generics here instead of two functions
void rebuild_kd_tree(bool is_initialized, bool copy_stuff, bool only_valid_in_curframe){
  //cout << "rebuild status " << status <<endl;
  if(is_initialized){
    delete local_priority_kdTree;
  }
  build_tree(LOCAL);
  
  if(copy_stuff){

    std::map<int, shared_ptr<vertex>> graph_points_loc = local_priority_kdTree_coordinates->get_vertices();
    //cout << "rebuild loc " << graph_points_loc.size() << endl;
    for (std::map<int, shared_ptr<vertex>>::iterator iterator = graph_points_loc.begin(); iterator != graph_points_loc.end(); iterator++){
      if(!only_valid_in_curframe){
        add_to_tree(iterator->second->get_location_coordinates(), iterator->second->get_index(), local_priority_kdTree);
      }
      else{
        if(!iterator->second->is_inactive()){
          add_to_tree(iterator->second->get_location_coordinates(), iterator->second->get_index(), local_priority_kdTree);
        }
      } 
    }
    //cout << "comparison: kd tree " << local_priority_kdTree->size << " data struct " << graph_points_loc.size() << endl;  
  }
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

bool test_path_noncolliding_static(shared_ptr<Path> path){

  shared_ptr<vertex> cur = path->get_beginning_node();
  shared_ptr<vertex> next = cur->get_child_pointer().lock();
  while(true){
    if(isnan(cur->get_location_coordinates()[0]) || isnan(cur->get_location_coordinates()[1]) || isnan(cur->get_location_coordinates()[2])){
      cout << "TESTING ERROR: NaN! " << endl;
      create_segfault(); 
    }
    if(is_in_obstacle_custom_frame(cur->get_location_coordinates(), cur->get_radius(), DONT_CHECK_WITH_BLOCKING_SPHERES, cur->get_first_frame())){
      cout << "TESTING ERROR: NODE IN PATH IS COLLIDING!" << endl;
      cout << "cur " << cur->get_index() << endl; 
      path->print_path();
      return false;
    }
    if(next->get_first_frame() != cur->get_last_valid_frame()){
      cout << "TESTING ERROR: PROBE IS NOT A FUCKING TARDIS!" << endl;
      cout << "cur " << cur->get_index() << endl; 
      path->print_path();
      return false;
    }

    if(next->get_index() == path->get_endpoint_index()){
      if(isnan(next->get_location_coordinates()[0]) || isnan(next->get_location_coordinates()[1]) || isnan(next->get_location_coordinates()[2])){
        cout << "TESTING ERROR: NaN! " << endl;
        create_segfault(); 
      }
      if(is_in_obstacle_custom_frame(next->get_location_coordinates(), next->get_radius(), DONT_CHECK_WITH_BLOCKING_SPHERES, next->get_first_frame())){
        cout << "TESTING ERROR: NODE IN PATH IS COLLIDING!" << endl;
        return false;
      }
      break;
    }

    cur = next;
    next = next->get_child_pointer().lock();
  } 

  return true;
} 
 
