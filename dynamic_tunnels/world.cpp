#include <iostream>
#include <fstream>
#include <sstream> 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <vector>

#include "world.h"
#include "ozcollide/ozcollide.h"
#include "MPNN/include/DNN/multiann.h"
#include "MPNN/include/DNN/ANN.h"

using namespace std;
using namespace ozcollide;
void create_spheres(const std::vector<Ball*> &balls);
vector<Ball*>& add_ball(double x, double y, double z, double radius, vector<Ball*> &balls);
vector<Ball*>& add_ball(double* loc_coord, double radius, vector<Ball*> &balls);
void delete_protein_structure();
AABBTreeSphere* build_protein_structure(const vector<Ball*> &balls);

AABBTreeSphere* protein_tree;
AABBTreeSphere* blocking_spheres_tree;
vector<Ball*> blocking_spheres;

bool TESTING_ENABLED = true;

double REPEATED_RUN_ITERATIONS_COEFFICIENT = 0.5;
int world_size_x = 100; int world_size_y = 100; int world_size_z = 100; int world_size_x_shift = 0; int world_size_y_shift = 0; int world_size_z_shift = 0;


//----parameters loaded from config file
int iterations;
double min_step;
double qix;             //starting position
double qiy;
double qiz;
string coords_file;     //file in which coordinates are located, option for static tunnel search, for dynamic search, individual frames should be placed in Frames directory and named [number of frame].txt
int number_of_frames;   
bool init_position_given; //program has the ability to find valid starting position by itself if not given
bool is_loaded_from_framesdir;
double probe_radius; //radius of sphere protein probe 
double test_sphere_radius; // radius of sphere used to find if point is inside protein structure
bool use_caver_dupcheck;
double inside_sampling_bias;
//---- variables for recording time spent in various methods
double nearest_neighbor_search_time = 0;
double tunnel_optimization_time = 0;
double centring_time = 0;
double collision_check_time = 0;
double duplication_check_time = 0;

//--- trees for nearest neighbor search
MPNN::MultiANN<double> *global_kdTree;
MPNN::MultiANN<double> *local_priority_kdTree;
Tree* local_priority_kdTree_coordinates; //exists only because dumbass MPNN doesn't support delete, kepts coordinates of local kd Tree
Tree* global_tree_points;                //to:do
std::vector <shared_ptr<Path>> paths;        //here are the fpund paths kept, all nodes in them are copied, therefore there isn't problem with deleting
                                  //their nodes everywhere. However, that doesn't protect you, if you want to access the node in main structure!
                                  //All nodes have related indices in main structure, however there is no guarantee that the ones in main structure
                                  //weren't deleted!

//--- identificators of global/local kd_tree for easier reading
const int NONE = 0;
const int GLOBAL = 1; 
const int LOCAL = 2;
const int BOTH = 3;
//--- current index serving as key inside kdtree
int tree_index = 0;
const int dimension = 3;


//---- global variables used in whole program
int cur_frame = 1; //frame which is currently used
bool rand_initialized = false; // for checking whether random number generator is initialized;

void init_random(){
  if(!rand_initialized){
    srand(time(NULL));
    rand_initialized = true;
  }
}

int get_current_frame(){
  return cur_frame;
}

int get_frames_count(){
  return number_of_frames;
}

void print_frame_borders(){
  cout << "printing frame borders" << endl;
  cout << "x: " << world_size_x_shift << " : " << world_size_x + world_size_x_shift << endl;
  cout << "y: " << world_size_y_shift << " : " << world_size_y + world_size_y_shift << endl;
  cout << "z: " << world_size_z_shift << " : " << world_size_z + world_size_z_shift << endl;
  cout << "end printing frame borders" << endl; 
}
AABBTreeSphere* get_blocking_spheres_tree(){
  //cout << "num l out" << blocking_spheres_tree->getNbLeafs() << endl;
  return blocking_spheres_tree;
}

bool is_init_position_given(){
  return init_position_given;
}
int get_blocking_spheres_tree_size(){
  return blocking_spheres_tree->getNbLeafs();
}
void print_blocking_spheres(){
  for(int i = 0; i < blocking_spheres.size(); i++){
    //cout << blocking_spheres[i]->location_coordinates[0] << "  " << blocking_spheres[i]->location_coordinates[1] << "  " << blocking_spheres[i]->location_coordinates[2] << endl;
  }
  //cout << endl;
}
int num_blocking_spheres(){
  return blocking_spheres.size();
}

void load_parameters(string config_filename){
  ifstream file(config_filename);
  string line;
  string value;
  string dummy;

  file >> dummy; file >> value;
  min_step = stod(value, NULL);

  file >> dummy; file >> value;
  iterations = stoi(value, NULL);

  file >> dummy; file >> value;
  probe_radius = stod(value, NULL);

  file >> dummy; file >> value;
  test_sphere_radius = stod(value, NULL);

  file >> dummy; file >> value;
  if(value == "is_not_given"){
    init_position_given = false;
  } else {
    qix = stod(value, NULL);
    init_position_given = true;
  } 

  file >> dummy; file >> value;
  if(init_position_given){
    qiy = stod(value, NULL);
  }

  file >> dummy; file >> value;
  if(init_position_given){
    qiz = stod(value, NULL);
  }

  file >> dummy; file >> value;
  coords_file = value;

  file >> dummy; file >> value;
  if(!(file.eof()) && stoi(value,NULL) > 0){
    number_of_frames = stoi(value, NULL);
    is_loaded_from_framesdir = true;
    cout << "frames: " << number_of_frames << endl; 
  } else {
    number_of_frames = 1;
    is_loaded_from_framesdir = false;
  }

  file >> dummy; file >> value;
  if(!(file.eof())){
    int dupcheck = stoi(value, NULL);
    if(dupcheck == 1){
      use_caver_dupcheck = true;
    } else {
      use_caver_dupcheck = false;
    }
  } else {
    use_caver_dupcheck = false;
  }
  file >> dummy; file >> value;
  if(!(file.eof())){
    inside_sampling_bias = stod(value, NULL);
    
  } else {
    inside_sampling_bias = 3;
  }


 
/*
  cout << min_step << endl;
  cout << iterations << endl;
  cout << probe_radius << endl;
  cout << test_sphere_radius << endl;
  cout << qix << endl;
  cout << qiy << endl;
  cout << qiz << endl;
*/
  file.close();

  
 
}
vector<Ball*> atoms_from_file(string coords_file){
  vector <Ball*> balls; 
  ifstream myfile(coords_file);

  string line;

  getline(myfile,line);
  double coordx; double coordy; double coordz; double radius;
  int counter = 0;
  double positive_border_X = -DBL_MAX; double negative_border_X = DBL_MAX;
  double positive_border_Y = -DBL_MAX; double negative_border_Y = DBL_MAX;
  double positive_border_Z = -DBL_MAX; double negative_border_Z = DBL_MAX;
  while(getline(myfile, line)){
    counter++; 
    myfile >> coordx;
    myfile >> coordy;
    myfile >> coordz;
    myfile >> radius;
    add_ball(coordx, coordy, coordz, radius, balls);

    if(coordx > positive_border_X){
      positive_border_X = coordx;
    }
    if(coordx < negative_border_X){
      negative_border_X = coordx;
    }
    if(coordy > positive_border_Y){
      positive_border_Y = coordy;
    }
    if(coordy < negative_border_Y){
      negative_border_Y = coordy;
    }
    if(coordz > positive_border_Z){
      positive_border_Z = coordz;
    }
    if(coordz < negative_border_Z){
      negative_border_Z = coordz;
    }
  }
  world_size_x = positive_border_X - negative_border_X + 4 * test_sphere_radius; world_size_x_shift = negative_border_X - 2 * test_sphere_radius;
  world_size_y = positive_border_Y - negative_border_Y + 4 * test_sphere_radius; world_size_y_shift = negative_border_Y - 2 * test_sphere_radius;
  world_size_z = positive_border_Z - negative_border_Z + 4 * test_sphere_radius; world_size_z_shift = negative_border_Z - 2 * test_sphere_radius;
  return balls;
}

void init_protein_struct(){
  //static tunnel detection problem assumed, uses file specified in config
  vector<Ball*> atoms;
  if(!is_loaded_from_framesdir){
    atoms = atoms_from_file(coords_file);
  } else {
    stringstream ss;
    ss << "Frames/" << cur_frame << ".txt";
    cout << ss.str() << endl;
    atoms = atoms_from_file(ss.str());
  }

  protein_tree = build_protein_structure(atoms);
  for(int i = 0; i < atoms.size(); i++){
    delete atoms[i];   
  }
  atoms.clear();
}

void run_next_frame(){
  cur_frame++;
  if(cur_frame == 2){  //expecting small inter-frame changes, it doesn't make much sense to run as many iterations as in the original
    iterations *= REPEATED_RUN_ITERATIONS_COEFFICIENT;
  }
  delete_protein_structure();
  init_protein_struct();
}


void rebuild_protein_structure(double* new_sphere_coords, bool add_sphere, string file_name){
  delete_protein_structure();
  vector<Ball*> atoms = atoms_from_file(file_name);
 // cout << new_sphere_coords[0] << " " << new_sphere_coords[1] << " " << new_sphere_coords[2] << endl;
  //if(add_sphere)add_ball(new_sphere_coords[0], new_sphere_coords[1], new_sphere_coords[2], test_sphere_radius, blocking_spheres);
  for(int i = 0; i < blocking_spheres.size(); i++){
    atoms.push_back(blocking_spheres[i]);
  }
  protein_tree = build_protein_structure(atoms);
}
void rebuild_blocking_spheres_structure(double* obstacle_loccoord, double radius){
  if(blocking_spheres.size() != 0) delete blocking_spheres_tree;
 
  Ball* obstacle = new Ball(obstacle_loccoord, radius);
  blocking_spheres.push_back(obstacle);
  
  blocking_spheres_tree = build_protein_structure(blocking_spheres);
  
  //mem error if there is only one sphere in tree for unknown reason
  if(blocking_spheres.size() == 1){
    rebuild_blocking_spheres_structure(obstacle_loccoord, radius);
  }
  //cout << "num l out" << blocking_spheres_tree->getNbLeafs() << endl;
}


vector<Ball*>& add_ball(double x, double y, double z, double radius, vector<Ball*> &balls){
  double* obstacle_loccoord = new double[3]; obstacle_loccoord[0] = x; obstacle_loccoord[1] = y; obstacle_loccoord[2] = z;
  Ball* obstacle = new Ball(obstacle_loccoord, radius);
  balls.push_back(obstacle);
  return balls;
}
vector<Ball*>& add_ball(double* loc_coord, double radius, vector<Ball*> &balls){
  Ball* obstacle = new Ball(loc_coord, radius);
  balls.push_back(obstacle);
  return balls;
}


AABBTreeSphere* build_protein_structure(const vector<Ball*> &balls){
  Sphere* spheres = new Sphere[balls.size()];

  for(int i=0;i<(int)balls.size(); i++) {
    spheres[i].center.x = balls[i]->location_coordinates[0];
    spheres[i].center.y = balls[i]->location_coordinates[1];
    spheres[i].center.z = balls[i]->location_coordinates[2];
    spheres[i].radius = balls[i]->radius;
  }
  
  AABBTreeSphere_Builder builder;
  AABBTreeSphere* tree = builder.build(balls.size(),spheres);
  delete [] spheres;
  spheres = NULL;
  return tree;
}

void delete_protein_structure(){
  delete protein_tree;
}

void delete_blocking_spheres(){
  //use shared ptr! it's impossible to hold on which one isn't already deleted
  blocking_spheres.clear();
  delete blocking_spheres_tree;
}