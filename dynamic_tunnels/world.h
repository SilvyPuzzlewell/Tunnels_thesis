#ifndef PARAM
#define PARAM
#include <stdbool.h>
#include <vector>
#include <map>
#include <memory>
#include <string>
#include "ozcollide/ozcollide.h"
#include "MPNN/include/DNN/multiann.h"
#include "MPNN/include/DNN/ANN.h"


#define M_PI 3.14159265358979323846

using namespace std;


extern int world_size_x; extern int world_size_y; extern int world_size_z; extern int world_size_x_shift; extern int world_size_y_shift; extern int world_size_z_shift;

extern bool TESTING_ENABLED;
//initial position and goal coordinates
extern double qix;
extern double qiy;
extern double qiz;

//algorithm parameters
extern double min_step; // maximum distance of new node from its parent
extern int iterations;
extern double min_distance_to_goal;
extern bool use_caver_dupcheck;
extern double inside_sampling_bias;
extern double MIN_VALID_INTERTUNNEL_DISTANCE;
extern bool RESETED_TREE_MODE;

//variables spend for recording time in various methods
extern double nearest_neighbor_search_time;
extern double tunnel_optimization_time;
extern double centring_time;
extern double collision_check_time;
extern double duplication_check_time;

//trees for nearest neighbor search
extern MPNN::MultiANN<double> *global_kdTree;
extern MPNN::MultiANN<double> *local_priority_kdTree;





class vertex {
  private:
  //vector<int> valid_frames; // frames in which this vertex is valid, should be sorted from lowest
  int exists_in_path;        // pointer into path in which the index; paths doesn't have explicit indexing, index into vector "paths" is used, this method expects paths in "paths" to not be deleted, if they are, this will cause crash
                            // or invalid results!
  vector<weak_ptr<vertex>> children_pointers;
  weak_ptr<vertex> parent_pointer;

  bool local;
  double* location_coordinates;
  int index;
  double radius;
  int frame_index; //used only for path vertices, denotes frame in which the vertex is in the path

  int first_frame;
  int last_valid_frame;
  bool inactive;

  public:
  static int vertex_counter;                 //predelat na class a dopsat destruktor
  
  vertex(double* location_coordinates_raw, shared_ptr<vertex> parent_pointer, int index, double radius, int first_frame, int last_valid_frame, bool local);
  vertex(double* location_coordinates_raw, int index, double radius, int first_frame, int last_valid_frame, bool local);
  ~vertex();

  //int find_previous_frame(int frame);
  //int get_last_valid_frame();
  shared_ptr<vertex> copy();
  shared_ptr<vertex> copy_without_structure_pointers();
  double* copy_coordinates();
  void set_coordinates(double* location_coordinates_raw);
  //void add_valid_frame(int frame);
  int is_in_path();
  bool is_in_some_valid_path();
  void set_in_path(int index);

  int get_index();
  void set_index(int index);

  void add_child_pointer(shared_ptr<vertex> child, bool print_output);
  bool is_child_local(int index);
  int get_children_count();
  void delete_children();

  bool is_local();
  bool is_parent_local();
  void make_global();
  int get_parent_index();
  weak_ptr<vertex> get_parent_pointer();
  void set_parent_pointer(shared_ptr<vertex> parent_pointer);

  int get_child_index(); //used for path nodes
  void set_child_pointer(shared_ptr<vertex> child_pointer);//used for path nodes
  weak_ptr<vertex> get_child_pointer();


  int get_child_index(int index); //for iterating over children_pointers vector, index is into the vector
  weak_ptr<vertex> get_child_pointer(int index);

  void delete_child(int index); //index is into the structure, not node index!

  //void set_frame_index(int frame_index);
  //int get_frame_index(); 
  void set_radius(double radius);
  double get_radius();
  double* get_location_coordinates();

  weak_ptr<vertex> get_child_pointer_null_permisive();
  weak_ptr<vertex> get_parent_pointer_null_permisive();

  //void set_first_frame(int first_frame);
  void set_last_valid_frame(int last_valid_frame);
  void set_inactive(bool inactive);
  int get_first_frame();
  int get_last_valid_frame();
  bool is_inactive();
  bool is_in_more_frames(); //shortcut for detecting time leaps, useful in tunnel optimization
  bool is_valid_in_frame(int frame);
  
}; 


class Tree {
 protected:
 std::map<int, shared_ptr<vertex>> path_vertices;
 public:
  
  ~Tree();
  shared_ptr<vertex> get_node_pointer(int index);
  double* get_node_coordinates(int index);
  void add_node(shared_ptr<vertex> node);
  int get_parent_index(int index);
  void print_nodes();
  std::map<int, shared_ptr<vertex>>& get_vertices();
  int get_size();
  int get_children_count(int index);
  void delete_node(int index);
  void reset();
  void print_vertices();
};

class Path: public Tree { 
 private:
  vector<int> valid_frames;
  //---for fast duplication checks
  
  vector<int> N_counts;
  int beginning_index;
  int endpoint_index;
  //---
 public:
  double bottleneck;
  double length;

  vector<double*> N_representation;
  shared_ptr<vertex> operator[](std::size_t idx); 

  Path(int beginning_index, int endpoint_index, int current_frame);
  ~Path();
  void add_valid_frame(int index);
  void erase_node_with_reconnecting(int index);
  void add_node(shared_ptr<vertex> node);
  int get_last_valid_frame();
  void add_N_count(int N_count);
  void add_N_point(double* N_point);
  bool check_N_representation_balance();
  double* get_N_point(int index);
  void reset_N_representation();
  void erase_from_index(int index);
  bool is_node_endpoint(shared_ptr<vertex> node);

  int get_child_index(int index);
  int get_endpoint_index();
  void set_endpoint_index(int endpoint_index);
  int get_beginning_index();
  void set_beginning_index(int beginning_index);

  void print_path();

  shared_ptr<vertex> get_endpoint_node();
  shared_ptr<vertex> get_beginning_node();
};

//todo - make world class

extern int get_current_frame();
extern void print_frame_borders();
extern bool is_init_position_given();
extern void init_random();
extern int num_blocking_spheres();
extern void init_protein_struct();
extern ozcollide::AABBTreeSphere* get_blocking_spheres_tree();
extern ozcollide::AABBTreeSphere* get_frame_protein_tree(int frame);
extern int get_blocking_spheres_tree_size();
extern void print_blocking_spheres();
extern void rebuild_protein_structure(double* new_sphere_coords, bool add_sphere);
extern void rebuild_blocking_spheres_structure(double* obstacle_loccoord, double radius);
extern void load_parameters(string config_filename);
extern int get_frames_count();
extern void run_next_frame();

extern void print_ownership(std::map<int, shared_ptr<vertex>> map, string message);

//for occollide
class Ball{
 public:
  Ball(double x, double y, double z, double radius);
  ~Ball();
  
  static int counter;
  double* location_coordinates;
  double radius;
};

extern double probe_radius;
extern double test_sphere_radius;

extern Tree* local_priority_kdTree_coordinates; //exists only because dumbass MPNN doesn't support delete, kepts coordinates of local kd Tree
extern Tree* global_tree_points;                //to:do
extern vector <shared_ptr<Path>> paths;        //here are the fpund paths kept, all nodes in them are copied, therefore there isn't problem with deleting
                                  //their nodes everywhere. However, that doesn't protect you, if you want to access the node in main structure!
                                  //All nodes have related indices in main structure, however there is no guarantee that the ones in main structure

//--- identificators of global/local kd_tree for easier reading
extern const int NONE;
extern const int GLOBAL; 
extern const int LOCAL;
extern const int BOTH;
//---

extern int tree_index;
extern const int dimension; //kd tree dimension
extern bool is_local_tree_initialized;



#endif