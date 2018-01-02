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



vertex::vertex(double* location_coordinates_raw, shared_ptr<vertex> parent_pointer, int index, double radius, int first_frame, int last_valid_frame, bool local){
    double* location_coordinates = new double[3]; location_coordinates[0] = location_coordinates_raw[0]; location_coordinates[1] = location_coordinates_raw[1]; location_coordinates[2] = location_coordinates_raw[2];
	this->location_coordinates = location_coordinates;
	this->parent_pointer = parent_pointer;
	this->radius = radius;
	this->index = index;
	this->local = local;
	vertex_counter++;
	if(first_frame != -1){ //-1 is used when for constructor in copying vertex, because we don't want to add new valid frame to it
		this->first_frame = first_frame;
	}
	if(last_valid_frame != -1){ //-1 is used when for constructor in copying vertex, because we don't want to add new valid frame to it
		this->last_valid_frame = last_valid_frame;
	}
	exists_in_path = -1; //newly created vertex doesn't exist in any path
	frame_index = 0; //initialization without path frame denotion
	inactive = false;
}

vertex::vertex(double* location_coordinates_raw, int index, double radius, int first_frame, int last_valid_frame, bool local){
    double* location_coordinates = new double[3]; location_coordinates[0] = location_coordinates_raw[0]; location_coordinates[1] = location_coordinates_raw[1]; location_coordinates[2] = location_coordinates_raw[2];
	this->location_coordinates = location_coordinates;
	//parent_pointer = NULL;
	this->radius = radius;
	this->index = index;
	this->local = local;
	vertex_counter++;
	if(first_frame != -1){ //-1 is used when for constructor in copying vertex, because we don't want to add new valid frame to it
		this->first_frame = first_frame;
	}
	if(last_valid_frame != -1){ //-1 is used when for constructor in copying vertex, because we don't want to add new valid frame to it
		this->last_valid_frame = last_valid_frame;
	}
	exists_in_path = -1; //newly created vertex doesn't exist in any path
	frame_index = 0; //initialization without path frame denotion
	inactive = false;
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

	return ret;
}

shared_ptr<vertex> vertex::copy_without_structure_pointers(){
	shared_ptr<vertex> ret = make_shared<vertex>(this->location_coordinates,this->index, this->radius, -1, -1, this->local);
	ret->first_frame = this->first_frame;
	ret->last_valid_frame = this->last_valid_frame;
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



int vertex::is_in_path(){
	return exists_in_path;
}


bool vertex::is_in_some_valid_path(){
	return this->exists_in_path == -1 ? false :  true;
}
void vertex::set_in_path(int new_index){
	this->exists_in_path = new_index;
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
	return children_pointers[index].lock()->is_local();
}
bool vertex::is_parent_local(){
	return parent_pointer.lock()->is_local();
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
	if(parent_pointer.expired()){
		return -1;
	}
	return parent_pointer.lock()->get_index();
}

int vertex::get_child_index(){
	if(children_pointers.size() != 1){
		cerr << "vertex::get_child_index():wrong use of vertex::get_child_index!" << endl;
		exit(0);
		return -1;
	}
	return (*children_pointers.begin()).lock()->get_index();
}

weak_ptr<vertex> vertex::get_child_pointer(int index){
	if(children_pointers.size() == 0){
		cerr << "vertex:: get_child_pointer error - no child!" << endl;
		//can't use null
		return weak_ptr<vertex>(); 
	}
	weak_ptr<vertex> ret = *(children_pointers.begin() + index);
	if(ret.expired()){
		cerr << "vertex:: get_child_pointer error - returning null pointer!" << endl;
	}
	return ret;
}

weak_ptr<vertex> vertex::get_child_pointer(){
	if(children_pointers.size() == 1){
		return children_pointers[0];
	} else if(children_pointers.size() > 1){
		cerr << "vertex::get_child_pointer(): children pointers are larger than one, illegal use" << endl;
		sleep(10);
		return children_pointers[0];
	} else {
		cerr << "vertex::get_child_pointer(): children pointers are empty" << endl;
		create_segfault();
		//return children_pointers[5];
		
	}

}

weak_ptr<vertex> vertex::get_child_pointer_null_permisive(){
	if(children_pointers.size() == 1){
		return children_pointers[0];
	} else if(children_pointers.size() > 1){
		cerr << "vertex::get_child_pointer(): children pointers are larger than one, illegal use" << endl;
		sleep(10);
		return children_pointers[0];
	}

	return weak_ptr<vertex>();	
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



int vertex::get_child_index(int index){
	if(children_pointers.size() == 0){
		cerr << "vertex::get_child_index(int index): attempt to get child of childless node!" << endl;
	}
	return children_pointers[index].lock()->get_index();
}

weak_ptr<vertex> vertex::get_parent_pointer(){
	if(parent_pointer.expired()){
		cerr << "vertex::get_parent_pointer(): program is trying to return null parent pointer" << endl;
	}
	return parent_pointer;
}

weak_ptr<vertex> vertex::get_parent_pointer_null_permisive(){
	return parent_pointer;
}



//void vertex::set_first_frame(int first_frame){
// 	this->first_frame = first_frame;
//}
void vertex::set_last_valid_frame(int last_valid_frame){
	this->last_valid_frame = last_valid_frame;
}
void vertex::set_inactive(bool inactive){
	this->inactive = inactive;
}
int vertex::get_first_frame(){
	return first_frame;
}
int vertex::get_last_valid_frame(){
	return last_valid_frame;
}

bool vertex::is_inactive(){
	return inactive;
} //shortcut for detecting time leaps, useful in tunnel optimization
bool vertex::is_in_more_frames(){
	return last_valid_frame - first_frame == 0 ? false : true;
}
bool vertex::is_valid_in_frame(int frame){
	//sanity check 
	if(first_frame > last_valid_frame){
		cout << "OBSTACLE.cpp, is_valid_in_frame insane " << endl;
		create_segfault();
	}
	if(frame <= last_valid_frame && frame >= first_frame){
		return true;
	} else {
		return false;
	}
} 


Path::Path(int beginning_index, int endpoint_index, int current_frame): beginning_index(beginning_index), endpoint_index(endpoint_index){
	valid_frames.push_back(current_frame);
}


void print_ownership(std::map<int, shared_ptr<vertex>> map, string message){
	int counter = 0;
	cout << endl;
	for(std::map<int, shared_ptr<vertex>>::iterator iterator = map.begin(); iterator != map.end(); iterator++){
	counter++;
	
	if(counter > 6){
		break;
	}
	
    cout << message << iterator->second.use_count() << endl;
   }
   cout<<endl;
}

Path::~Path(){

	reset_N_representation();
}


void Path::add_valid_frame(int frame){
	valid_frames.push_back(frame);
}

bool Path::is_node_endpoint(shared_ptr<vertex> node){
	return node->get_index() == this->get_endpoint_index();
}

shared_ptr<vertex> Path::get_beginning_node(){
	return path_vertices[this->get_beginning_index()];
}

shared_ptr<vertex> Path::get_endpoint_node(){
	return path_vertices[this->get_endpoint_index()];
}

shared_ptr<vertex> Path::operator[](std::size_t idx){
	return path_vertices[idx];
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

void Tree::print_vertices(){
	print_map_indices(path_vertices, "");
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
	weak_ptr<vertex> parent_pointer = path_vertices[index]->get_parent_pointer_null_permisive();
	weak_ptr<vertex> child_pointer = path_vertices[index]->get_child_pointer_null_permisive();
	if(!parent_pointer.expired()){
	//	cout << "debug parent_index " << parent_pointer->get_index() << endl;
	//	cout << "debug parent pointer " << parent_pointer << " local " << parent_pointer->is_local() <<endl;

	}
	if(!child_pointer.expired()){
	//	cout << "debug child_index " << child_pointer->get_index() << endl;
	//	cout << "debug child pointer " << child_pointer << " local " <<child_pointer->is_local() <<endl;
	}

    if(!parent_pointer.expired() && !child_pointer.expired()){
    //	cout << "debug = opt 1" << endl;
		parent_pointer.lock()->set_child_pointer(child_pointer.lock()); // connect parent node to it's child
		child_pointer.lock()->set_parent_pointer(parent_pointer.lock());
	} else if(parent_pointer.expired() && !child_pointer.expired()){
	//	cout << "debug = opt 2" << endl;
		cerr << "warning -- deleted root path node" << endl;
		child_pointer.lock()->set_parent_pointer(NULL);
		beginning_index = child_pointer.lock()->get_index();
	} else if(!parent_pointer.expired() && child_pointer.expired()){
	//	cout << "debug = opt 3" << endl;
	//	cout << "parent pointer " << parent_pointer << endl;
		parent_pointer.lock()->set_child_pointer(NULL);
		endpoint_index = parent_pointer.lock()->get_index();
	} else {
	//	cout << "debug = opt 4" << endl;
		cerr << "Path::erase_node_with_reconnecting(int index): deleting from one atom path!!" << endl;
		exit(0);
	}	
	path_vertices.erase(index);
}


void Path::add_N_point(double* N_point){ 	
	N_representation.push_back(copy_vector(N_point));	
}


void Path::reset_N_representation(){
	delete_vector(this->N_representation);
	N_representation.clear();
}

double* Path::get_N_point(int index){
	if(index >= N_representation.size()){
		cerr << "fatal error, asked for higher N in checking duplications" << endl;
		create_segfault();
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

void Path::print_path(){
	shared_ptr<vertex> cur = this->get_beginning_node();
	shared_ptr<vertex> next = cur->get_child_pointer().lock();
	while(true){
		cout << "index " << cur->get_index() << endl;
		cout << "first frame " << cur->get_first_frame() << endl;
		cout << "last valid frame " << cur->get_last_valid_frame() << endl;
		cout << "is in path " << cur->is_in_path() << endl;
		cout << "coordinates" << endl;
		print_vector(cur->get_location_coordinates());
		cout << "distance to next " << compute_metric_eucleidean(cur->get_location_coordinates(), next->get_location_coordinates()) << endl;

		if(next->get_index() == endpoint_index){
			cout << "endpoint " << endl;
			cout << "index " << next->get_index() << endl;
			cout << "first frame " << next->get_first_frame() << endl;
			cout << "last valid frame " << next->get_last_valid_frame() << endl;
			cout << "is in path " << next->is_in_path() << endl;
			cout << "coordinates" << endl;
			print_vector(next->get_location_coordinates());
			break;
		}

		cur = next;
		next = next->get_child_pointer().lock();
	} 
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
	return path_vertices[index]->get_child_pointer(0).lock()->get_index();
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

