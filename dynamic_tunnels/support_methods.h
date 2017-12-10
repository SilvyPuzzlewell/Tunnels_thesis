
#ifndef SUPPORT_METHODS
#define SUPPORT_METHODS

#include <vector>
#include <string>
#include <memory>
#include "world.h"

const int ADDITION = 0;
const int SUBTRACTION = 1;
double* add_vectors(double* first, double* second, int mode);
double* copy_vector(double* vector);
double dot_product(double* first, double* second);
void multiply_vector(double* array, double constant);
double vector_length(double* vector);
void normalize_vector(double* vector,double normalised_length);
double compute_metric_eucleidean(double* coord1, double* coord2, int dimension);
double compute_metric_eucleidean(double* coord1, double* coord2);
void print_vector(double* vector);
void divide_vector(double* vector, double constant);
int return_min(int i1, int i2);
void delete_vertex_pointer_vector(std::vector<shared_ptr<vertex>>& vector);
void delete_vertex_pointer_vector(std::map<int, shared_ptr<vertex>>& vector);
void print_map_indices(std::map<int ,shared_ptr<vertex>>& map, std::string message);
int signum(double value);
bool contains_key(int index, std::map<int, shared_ptr<vertex>>& map);
double compute_metric_eucleidean(array<double, 3> coord1, array<double, 3> coord2);
array<double, 3> copy_coordinate_pointer_to_array(double* pointer);
double* copy_array_to_coordinate_pointer(array<double, 3> array);
bool double_equals(double a, double b);


#endif