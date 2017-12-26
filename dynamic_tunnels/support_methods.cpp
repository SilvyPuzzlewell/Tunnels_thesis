#include <cmath>
#include <iostream>
#include <string>
#include <memory>
#include <signal.h>

#include "support_methods.h"


double* add_vectors(double* first, double* second, int mode){
  if(isnan(first[0])){
    create_segfault();
  }
  if(isnan(second[0])){
    create_segfault();
  }
  double* ret = new double[3]; //delet this later
  if(mode == ADDITION){
    ret[0] = first[0] + second[0]; ret[1] = first[1] + second[1]; ret[2] = first[2] + second[2];
  } else if(mode == SUBTRACTION){
    ret[0] = first[0] - second[0]; ret[1] = first[1] - second[1]; ret[2] = first[2] - second[2];
  }
  if(isnan(first[0])){
    create_segfault();
  }
  return ret;
}
int return_min(int i1, int i2){
  return i1 < i2 ? i1 : i2;
}
double* copy_vector(double* vector){
  double* ret = new double[3]; //delet this later
  double x_coord = vector[0]; double y_coord = vector[1]; double z_coord = vector[2];
  ret[0] = x_coord; ret[1] = y_coord; ret[2] = z_coord;
  if(isnan(ret[0])){
    create_segfault();
  }
  return ret;
}
//for sanity checks
double dot_product(double* first, double* second){
  return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}
//works only for arrays initialized by new with size 3!!
void multiply_vector(double* array, double constant){
  for(int i = 0; i < 3; i++){
    array[i] = array[i] * constant;
  }
  if(isnan(array[0])){
    create_segfault();
  }
}

double vector_length(double* vector){
  double distance_squared = 0;
  for(int i = 0; i < 3; i++){
    distance_squared += pow(vector[i], 2);
  }
  return sqrt(distance_squared);
}

//allocated with new
void normalize_vector(double* vector,double normalised_length){
    double length = vector_length(vector);
    if(length == 0){
      create_segfault();
    }
    if(isnan(vector[0])){
      create_segfault();
    }
    for(int i = 0; i < 3; i++){
      vector[i] = (vector[i] / length) * normalised_length;
    }
    if(isnan(vector[0])){
     create_segfault();
    }
}

void divide_vector(double* vector, double constant){
  for(int i = 0; i < 3; i++){
    vector[i] = vector[i] / constant;
  }
}

double compute_metric_eucleidean(double* coord1, double* coord2, int dimension){
  double distance_squared = 0;
  for(int i = 0; i < dimension; i++){
    distance_squared += pow(coord1[i]  - coord2[i], 2);
  }
  return sqrt(distance_squared);
}
//default in 3d space
double compute_metric_eucleidean(double* coord1, double* coord2){
  double distance_squared = 0;
  for(int i = 0; i < 3; i++){
    distance_squared += pow(coord1[i]  - coord2[i], 2);
  }
  return sqrt(distance_squared);
}

double compute_metric_eucleidean(array<double, 3> coord1, array<double, 3> coord2){
  double distance_squared = 0;
  for(int i = 0; i < 3; i++){
    distance_squared += pow(coord1[i]  - coord2[i], 2);
  }
  return sqrt(distance_squared);
}


void print_vector(double* vector){
  std::cout << " printing vector "<< vector[0] << " " << vector[1] << " " << vector[2] <<std::endl;
}

void delete_vertex_pointer_vector(std::vector<shared_ptr<vertex>>& vector){
  while(vector.size() != 0){
    vector.pop_back();
  } 
}

void delete_vertex_pointer_vector(std::map<int, shared_ptr<vertex>>& vector){
  for(std::map<int, shared_ptr<vertex>>::iterator iterator = vector.begin(); iterator != vector.end(); iterator++){
   vector.clear();
  }
}
void print_map_indices(std::map<int ,shared_ptr<vertex>>& map, std::string message){
  std::cout << std::endl << "printing map coordinates " << std::endl << message << std::endl;
  for(std::map<int, shared_ptr<vertex>>::iterator iterator = map.begin(); iterator != map.end(); iterator++){
    std::cout << " " << iterator->first;
   }
   std::cout << std::endl;
}

int signum(double value){
  return value < 0 ? -1 : 1;
}

bool contains_key(int index, std::map<int, shared_ptr<vertex>>& map){
  std::map<int, shared_ptr<vertex>>::iterator it;
  it = map.find(index);
  return it != map.end();
}

//stupid hack, I discovered the array thing later and i'm lazy ya know, I don't want to redo whole code now
array<double, 3> copy_coordinate_pointer_to_array(double* pointer){
  return array<double, 3>{pointer[0], pointer[1], pointer[2]};  
}

//as before, AND DON'T FORGET THE MEMORY MANAGEMENT
double* copy_array_to_coordinate_pointer(array<double, 3> array){
  double* coordinates = new double[3]; 
  coordinates[0] = array[0]; coordinates[1] = array[1]; coordinates[2] = array[2];
  return coordinates;
}

double* create_direction_vector(shared_ptr<vertex> base, shared_ptr<vertex> target, double length){
  double* direction_vector = add_vectors(target->get_location_coordinates(), base->get_location_coordinates(), SUBTRACTION);  
  normalize_vector(direction_vector, length);
  return direction_vector;
}

bool double_equals(double a, double b)
{
    double epsilon = 0.00001;
    return std::abs(a - b) < epsilon;
}

//This is actually useful in debugging!
void create_segfault(){
  raise(SIGSEGV);
}

void delete_vector(std::vector<double*> v){
  //cout << "size pre " << v.size() << endl;
  while(v.size() != 0){
      double* deleted_ptr = v.back();
      delete [] deleted_ptr;
      v.pop_back();
    }
  //  cout << "size post " << v.size() << endl;
} 