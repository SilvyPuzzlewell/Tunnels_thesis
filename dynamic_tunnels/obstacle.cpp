#include <cmath>
#include "world.h"
#include "ozcollide/ozcollide.h"
#include "support_methods.h"
#include <iostream>
#include <string.h>
#include <memory>
#include <map>
#include <sys/stat.h>
#include <algorithm>
#include <utility> 
#include <unistd.h>
#include <string>
using namespace std;



vertex::vertex(double* location_coordinates_raw, shared_ptr<vertex> parent_pointer, int index, double radius, int cur_frame, bool local){
    double* location_coordinates = new double[3]; location_coordinates[0] = location_coordinates_raw[0]; location_coordinates[1] = location_coordinates_raw[1]; location_coordinates[2] = location_coordinates_raw[2];
	this->location_coordinates = location_coordinates;
	this->parent_pointer = parent_pointer;
	this->radius = radius;
	this->index = index;
	this->local = local;
	vertex_counter++;
	if(cur_frame != -1){ //-1 is used when for constructor in copying vertex, because we don't want to add new valid frame to it
		valid_frames.push_back(cur_frame);
	}
	exists_in_path = -1; //newly created vertex doesn't exist in any path
	frame_index = 0; //initialization without path frame denotion
}

vertex::vertex(double* location_coordinates_raw, int index, double radius, int cur_frame, bool local){
    double* location_coordinates = new double[3]; location_coordinates[0] = location_coordinates_raw[0]; location_coordinates[1] = location_coordinates_raw[1]; location_coordinates[2] = location_coordinates_raw[2];
	this->location_coordinates = location_coordinates;
	parent_pointer = NULL;
	this->radius = radius;
	this->index = index;
	this->local = local;
	vertex_counter++;
	if(cur_frame != -1){ //-1 is used when for constructor in copying vertex, because we don't want to add new valid frame to it
		valid_frames.push_back(cur_frame);
	}
	exists_in_path = -1; //newly created vertex doesn't exist in any path
	frame_index = 0; //initialization without path frame denotion
}

vertex::~vertex(){
	//std::cout << location_coordinates[0] << endl;
	delete [] location_coordinates;
	children_pointers.clear();
	vertex_counter--;
}
	shared_ptr<vertex> vertex::copy(){
	shared_ptr<vertex> ret = this->copy_without_structure_pointers();
	ret->children_pointers = this->children_pointers;
	ret->valid_frames = this->valid_frames;
	return ret;
}

shared_ptr<vertex> vertex::copy_without_structure_pointers(){
	shared_ptr<vertex> ret = make_shared<vertex>(this->location_coordinates,this->index, this->radius, -1, this->local);
	ret->valid_frames = this->valid_frames;
	//cout << "SHARED COPIED " << ret.use_count() << endl;
	return ret;
}
double* vertex::copy_coordinates(){
	double* ret = new double[3];
	double x_coord = this->location_coordinates[0]; double y_coord = this->location_coordinates[1]; double z_coord = this->location_coordinates[2];
	ret[0] = x_coord; ret[1] = y_coord; ret[2] = z_coord;
	return ret;
}

void vertex::set_coordinates(double* location_coordinates_raw){
	delete [] location_coordinates;
	double* location_coordinates = new double[3]; location_coordinates[0] = location_coordinates_raw[0]; location_coordinates[1] = location_coordinates_raw[1]; location_coordinates[2] = location_coordinates_raw[2];
	this->location_coordinates = location_coordinates;
	
}

void vertex::add_valid_frame(int frame){
	valid_frames.push_back(frame);
}

//expects sorted array, so it can find the largest passable frame
//it should return the next after the one passed in argument, therefore it should run in 1 time and should not require sophisticated search method
int vertex::find_previous_frame(int frame){
	for(int i = valid_frames.size() - 1; i >= 0; i--){
		if(valid_frames[i] < frame){
			return valid_frames[i];
		}
	}
	//cout << "fail " << valid_frames.size() << " " << frame << " " << valid_frames[0] << endl;
	return -1; //no previous frame exists
}
int vertex::get_last_valid_frame(){
	return valid_frames.back();
}

int vertex::is_in_path(){
	return exists_in_path;
}

void vertex::set_in_path(int new_index){
	this->exists_in_path = new_index;
}

bool vertex::is_valid_in_frame(int frame){
	//cout << "i" << endl;
	return std::binary_search(valid_frames.begin(), valid_frames.end(), frame);
}
int vertex::get_index(){
	return index;
}

void vertex::set_index(int index){
	this->index = index;
}

void vertex::add_child_pointer(shared_ptr<vertex> child, bool print_output){
	if(child == NULL){
		cerr << "vertex::add_child_pointer - adding null pointer!" << endl;
	}
	if(print_output){
		cout << "vertex::add_child_pointer adding child " << child->get_index() << " to " << index << endl;
	}
	children_pointers.push_back(child);
}

bool vertex::is_child_local(int index){
	return children_pointers[index]->is_local();
}
bool vertex::is_parent_local(){
	return parent_pointer->is_local();
}
void vertex::make_global(){
	local = false;
}
void vertex::delete_child(int index){
	children_pointers.erase(children_pointers.begin() + index);
 }

int vertex::get_children_count(){
	return children_pointers.size();
}


void vertex::delete_children(){
	children_pointers.clear();
}

bool vertex::is_local(){
	return local;
}

int vertex::get_parent_index(){
	if(parent_pointer == NULL){
		return -1;
	}
	return parent_pointer->get_index();
}

int vertex::get_child_index(){
	if(children_pointers.size() != 1){
		cerr << "vertex::get_child_index():wrong use of vertex::get_child_index!" << endl;
		exit(0);
		return -1;
	}
	return (*children_pointers.begin())->get_index();
}

shared_ptr<vertex> vertex::get_child_pointer(int index){
	if(children_pointers.size() == 0){
		cerr << "vertex:: get_child_pointer error - no child!" << endl;
		return NULL;
	}
	shared_ptr<vertex> ret = *(children_pointers.begin() + index);
	if(ret == NULL){
		cerr << "vertex:: get_child_pointer error - returning null pointer!" << endl;
	}
	return ret;
}

shared_ptr<vertex> vertex::get_child_pointer(){
	if(children_pointers.size() == 1){
		return children_pointers[0];
	} else if(children_pointers.size() > 1){
		cerr << "vertex::get_child_pointer(): children pointers are larger than one, illegal use" << endl;
		sleep(10);
		return children_pointers[0];
	} else {
		cerr << "vertex::get_child_pointer(): children pointers are empty" << endl;
		exit(0);
	}

}

shared_ptr<vertex> vertex::get_child_pointer_null_permisive(){
	if(children_pointers.size() == 1){
		return children_pointers[0];
	} else if(children_pointers.size() > 1){
		cerr << "vertex::get_child_pointer(): children pointers are larger than one, illegal use" << endl;
		sleep(10);
		return children_pointers[0];
	}

	return NULL;	
}


double vertex::get_radius(){
	return radius;
}

void vertex::set_radius(double radius){
	this->radius = radius;
}

double* vertex::get_location_coordinates(){
	return location_coordinates;
}

void vertex::set_parent_pointer(shared_ptr<vertex> parent_pointer){
	this->parent_pointer = parent_pointer;
}

void vertex::set_child_pointer(shared_ptr<vertex> child_pointer){
	if(children_pointers.size() > 1){
		cerr << "fail at set child in path, more than one pointer already" << endl;
	}
	children_pointers.clear();
	if(child_pointer != NULL){
		children_pointers.push_back(child_pointer);
	}
}

void vertex::set_frame_index(int frame_index){
	this->frame_index = frame_index;
}

int vertex::get_frame_index(){
	return frame_index;
}

int vertex::get_child_index(int index){
	if(children_pointers.size() == 0){
		cerr << "vertex::get_child_index(int index): attempt to get child of childless node!" << endl;
	}
	return children_pointers[index]->get_index();
}

shared_ptr<vertex> vertex::get_parent_pointer(){
	if(parent_pointer == NULL){
		cerr << "vertex::get_parent_pointer(): program is trying to return null parent pointer" << endl;
	}
	return parent_pointer;
}

shared_ptr<vertex> vertex::get_parent_pointer_null_permisive(){
	return parent_pointer;
}

Path::Path(int beginning_index, int endpoint_index, int current_frame): beginning_index(beginning_index), endpoint_index(endpoint_index){
	valid_frames.push_back(current_frame);
}


void print_ownership(std::map<int, shared_ptr<vertex>> map, string message){
	int counter = 0;
	for(std::map<int, shared_ptr<vertex>>::iterator iterator = map.begin(); iterator != map.end(); iterator++){
	counter++;
	if(counter > 6){
		break;
	}
    cout << message << iterator->second.use_count() << endl;
   }
}

Path::~Path(){
	//cout << "DELET PATH" << endl;

	/*
	map<int, shared_ptr<vertex>> tempVector;
	path_vertices.swap(tempVector);
	*/

  //---cyclic ownership shared_ptr memory leak happens without this 	
  for(std::map<int, shared_ptr<vertex>>::iterator iterator = path_vertices.begin(); iterator != path_vertices.end(); iterator++){
  	iterator->second->get_parent_pointer().reset();
  	iterator->second->get_child_pointer().reset();
   }
  //---


   //print_ownership(path_vertices, "BEFORE DELET");
	
	
	while(N_representation.size() != 0){
		double* deleted_ptr = N_representation.back();
		delete [] deleted_ptr;
		N_representation.pop_back();
	}
	
}


void Path::add_valid_frame(int frame){
	valid_frames.push_back(frame);
}
//it ia expected that nodes are copied to path, therefore deleting them does not cause mem corruption anywhere else
Tree::~Tree(){
   path_vertices.clear();
}

std::map<int, shared_ptr<vertex>>& Tree::get_vertices(){
	return path_vertices;
}

int Tree::get_size(){
	return path_vertices.size();
}


int Tree::get_children_count(int index){
	return path_vertices[index]->get_children_count();
}

void Tree::add_node(shared_ptr<vertex> node){
	path_vertices.insert(std::pair<int,shared_ptr<vertex>>(node->get_index(), node));
}

void Tree::delete_node(int index){
	//std::cout << "Tree: deleted node " << index << endl; 
	path_vertices.erase(index);
}

void Tree::reset(){
	path_vertices.clear();
}

void Path::add_node(shared_ptr<vertex> node){
	if(node->get_children_count() > 1){
		std::cerr << "fatal error - path node has more than one child";
		exit(0);
	}
	path_vertices.insert(std::pair<int,shared_ptr<vertex>>(node->get_index(), node));
}

int Path::get_last_valid_frame(){
	return valid_frames.back();
}

shared_ptr<vertex> Tree::get_node_pointer(int index){
	return this->path_vertices[index];
}


//parentless node is denoted by -1
int Tree::get_parent_index(int index){
	return this->path_vertices[index]->get_parent_index();
}
double* Tree::get_node_coordinates(int index){
	return this->path_vertices[index]->get_location_coordinates();
}

//wrapper function for deleting nodes from path without damaging it's structure
void Path::erase_node_with_reconnecting(int index){
	shared_ptr<vertex> parent_pointer = path_vertices[index]->get_parent_pointer_null_permisive();
	shared_ptr<vertex> child_pointer = path_vertices[index]->get_child_pointer_null_permisive();
	if(parent_pointer != NULL){
	//	cout << "debug parent_index " << parent_pointer->get_index() << endl;
	//	cout << "debug parent pointer " << parent_pointer << " local " << parent_pointer->is_local() <<endl;

	}
	if(child_pointer != NULL){
	//	cout << "debug child_index " << child_pointer->get_index() << endl;
	//	cout << "debug child pointer " << child_pointer << " local " <<child_pointer->is_local() <<endl;
	}

    if(parent_pointer != NULL && child_pointer != NULL){
    //	cout << "debug = opt 1" << endl;
		parent_pointer->set_child_pointer(child_pointer); // connect parent node to it's child
		child_pointer->set_parent_pointer(parent_pointer);
	} else if(parent_pointer == NULL && child_pointer != NULL){
	//	cout << "debug = opt 2" << endl;
		cerr << "warning -- deleted root path node" << endl;
		child_pointer->set_parent_pointer(NULL);
		beginning_index = child_pointer->get_index();
	} else if(parent_pointer != NULL && child_pointer == NULL){
	//	cout << "debug = opt 3" << endl;
	//	cout << "parent pointer " << parent_pointer << endl;
		parent_pointer->set_child_pointer(NULL);
		endpoint_index = parent_pointer->get_index();
	} else {
	//	cout << "debug = opt 4" << endl;
		cerr << "Path::erase_node_with_reconnecting(int index): deleting from one atom path!!" << endl;
		exit(0);
	}	
	path_vertices.erase(index);
}


void Path::add_N_point(double* N_point){ //added paths are never supposed to be deleted and N_points are limited only to them, hence no need for copying	
	N_representation.push_back(N_point);	
}


void Path::reset_N_representation(){
	N_representation.clear();
}

double* Path::get_N_point(int index){
	if(index >= N_representation.size()){
		cerr << "fatal error, asked for higher N in checking duplications" << endl;
		exit(0);
	} 
	else {
		return N_representation[index];
	}
}

int erase_index = 0;
void Path::erase_from_index(int index){
	erase_index++;
	//cout <<"debug = erase_index " << erase_index << endl;
	while(endpoint_index != index){
	//	this->print_nodes();
	//	cout << "debug " << endpoint_index << endl;
	//	cout <<"debug = erase pointer " << path_vertices[endpoint_index] << " local " <<path_vertices[endpoint_index]->is_local()<< endl;
		erase_node_with_reconnecting(endpoint_index);
		if(endpoint_index == beginning_index){
			cout << "warning - erasing whole path! Possible that index in erase_from_index doesn't exist in path" << endl;
			break;
		}
	}
	erase_node_with_reconnecting(index);
}

int Path::get_endpoint_index(){
	return endpoint_index;
}

void Path::set_endpoint_index(int endpoint_index){
	this->endpoint_index = endpoint_index;
}
int Path::get_beginning_index(){
  	return beginning_index;
}
void Path::set_beginning_index(int beginning_index){
	this->beginning_index = beginning_index;
}
int Path::get_child_index(int index){
	return path_vertices[index]->get_child_pointer(0)->get_index();
}

void Tree::print_nodes(){
	cout << "--printing path--" << endl;
	for(std::map<int, shared_ptr<vertex>>::iterator iterator = path_vertices.begin(); iterator != path_vertices.end(); iterator++){
		cout << iterator->second->get_index() << " "; 
	}
	cout << "--printing path--end--" << endl;
}

Ball::Ball(double x, double y, double z, double radius){
	double* location_coordinates = new double[3]; location_coordinates[0] = x; location_coordinates[1] = y; location_coordinates[2] = z;
	this->location_coordinates = location_coordinates;
	this->radius = radius;
	counter++;
}

Ball::~Ball(){
	counter--;
	delete [] location_coordinates;
}

