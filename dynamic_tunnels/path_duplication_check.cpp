#include <cmath>
#include <iostream>
#include <climits>
#include <map>
#include <vector>
#include <memory>

#include "support_methods.h"
#include "world.h"



int N = 50; //default number of intervals
double INTERVAL_RATIO_CONSTANT = 15; //into how much intervals are the paths divided, based on maximum number of tunnel atoms in current run, example for shortest path with 1000 atoms, run will be divided into 100 intervalsw with value 10 
double Z_start_radius = 1;
double Z_end_radius = 1;



bool DEBUG_ENABLED = false;


//decides if node lies in the z interval sphere
bool check_if_in_sphere(double* coordinate ,double sphere_radius,double* sphere_coordinate){
	double distance = compute_metric_eucleidean(coordinate, sphere_coordinate, 3);
	//std::cout << "distance " << distance << std::endl; 
	return distance < sphere_radius ? true : false;
}
	
int find_extreme_point(std::map<int, std::shared_ptr<vertex>> tunnel, double* center){
	double max_distance = -1;
	int extreme_point = -1; //index into tunnel denoting the extreme point
	for(std::map<int, std::shared_ptr<vertex>>::iterator iterator = tunnel.begin(); iterator != tunnel.end(); iterator++){
		double cur_distance = compute_metric_eucleidean(iterator->second->get_location_coordinates(), center);

		if(cur_distance > max_distance){
			cur_distance = max_distance;
			extreme_point = iterator->first;
		} 
	}
	if(extreme_point == -1){ //this shoudn't ever happen
		std::cerr << "fatal error in path duplication check - extreme point not found!" << std::endl;
		exit(0);
	}
	return extreme_point;
}


//finds the largest distance from center to tunnel node, beginning at start_index and ending at end_index, other nodes are ignored, returns extreme points index in tunnel
//this one is used
int find_extreme_point(shared_ptr<Path> tunnel, double* center, int start_index, int end_index){
	double max_distance = -1;
	int extreme_point = -1; //index into tunnel denoting the extreme point
	std::map<int, std::shared_ptr<vertex>>& tunnel_map = tunnel->get_vertices();

	int index = start_index;

	int counter = 0;
	while(true){
		double cur_distance = compute_metric_eucleidean(tunnel_map[index]->get_location_coordinates(), center);
		if(DEBUG_ENABLED){
			cout << "FINDING EXTREME POINT: current distance " << cur_distance << " current index " << index << endl;
		}
		if(cur_distance > max_distance){
			max_distance = cur_distance;
			extreme_point = index;
		}
		
		if(index == end_index){
			break;
		}
		index = tunnel->get_child_index(index);
		//std::cout << "find_extreme_point " << cur_distance << std::endl; 
		counter++;
	}
	
	if(extreme_point == -1){ //this shoudn't ever happen
		std::cerr << "fatal error in path duplication check - extreme point not found!" << std::endl;
		exit(0);
	}
	return extreme_point;
}
	
//Uses Z_interval parameters for not counting parts around beginning or ending. If used, delegation into intervals will start at first unremoved node and end at last node before the removed ones
//at the end.
bool create_N_representation(shared_ptr<Path> tunnel){
	int cutted_part_end_index = tunnel->get_endpoint_index();
	int cutted_part_start_index = tunnel->get_beginning_index();
	std::map<int, std::shared_ptr<vertex>> tunnel_map = tunnel->get_vertices();
	double* tunnel_end_coordinate = tunnel->get_endpoint_node()->get_location_coordinates();
	double* tunnel_root_coordinate = tunnel->get_beginning_node()->get_location_coordinates();
	
    //these take into account from where the representation is computed, ie where the the tunnel starts after getting rid of Z_start interval
	int start_offset = 0;
	int end_offset = 0;

	std::shared_ptr<vertex> cur = tunnel->get_beginning_node();
	bool are_endings_cut = true;

	

	//gets rid of nodes in the z radii
	double tunnel_offset = DBL_MAX;
	int counter = 0;

  	while(true){
  		double cur_distance_from_root = compute_metric_eucleidean(cur->get_location_coordinates(), tunnel_root_coordinate);
  		if(cur_distance_from_root<= Z_start_radius) {
  			cutted_part_start_index = cur->get_index();
  			
  		} else if(cur_distance_from_root > Z_start_radius){
  			if(cur_distance_from_root < tunnel_offset){
  				tunnel_offset = cur_distance_from_root;
  			}
  		}
  		if(check_if_in_sphere(cur->get_location_coordinates(), Z_end_radius, tunnel_end_coordinate)){
  			cutted_part_end_index = cur->get_index();
  			break;
  		}
  		if(tunnel->is_node_endpoint(cur)){
  			break;
  		}
  		cur = cur->get_child_pointer().lock();
  			
	}

	//we need to use the index next to the last one colliding with z_interval sphere; 
	cutted_part_start_index = tunnel->get_node_pointer(cutted_part_start_index)->get_child_pointer().lock()->get_index();
	cutted_part_end_index = tunnel->get_node_pointer(cutted_part_end_index)->get_parent_pointer().lock()->get_index();

	if(cutted_part_start_index >= cutted_part_end_index){
		are_endings_cut = false;
		cutted_part_start_index = tunnel->get_beginning_index();
		cutted_part_end_index = tunnel->get_endpoint_index();
		tunnel_offset = 0;
	} 
	

	//---gets rid of nodes in the z radii

	//locates node lying farthest from the tunnel centre, this is used as basis for the size of the intervals
	int extreme_point = find_extreme_point(tunnel, tunnel_root_coordinate, cutted_part_start_index, cutted_part_end_index);
	double extreme_point_distance = compute_metric_eucleidean(tunnel_map[extreme_point]->get_location_coordinates(), tunnel_root_coordinate) - tunnel_offset;
	//coordinate of the first node used for interval computation
	double* tunnel_beginning_coordinate = tunnel_map[cutted_part_start_index]->get_location_coordinates();
	//we must subtract the distance from the tunnel's first to node to first node after getting rid of nodes within z_interval, so that we can start at zero, which is required
	//for dividing into intervals correctly
	//cout << "entering test" << endl;
	//print_vector(tunnel_map[tunnel->get_beginning_index()]->get_location_coordinates());

	int index = cutted_part_start_index;
	std::vector<int> N_counts(N, 0);
	//just init zero values
	std::vector<double*> N_representation;
	for(int i = 0; i < N; i++){
		double* coordinate = new double[3];
		N_representation.push_back(coordinate);
	}

	counter = 0;
	while(true){
		//test
		if(TESTING_ENABLED){
			if(compute_metric_eucleidean(tunnel_map[index]->get_location_coordinates(), tunnel_root_coordinate)  <= Z_start_radius && are_endings_cut){
				cout << "TEST ERROR: node in z_radius in duplication check" << endl;
				cout << "distance from tunnel root: " <<  compute_metric_eucleidean(tunnel_map[index]->get_location_coordinates(), tunnel_root_coordinate) << endl;
				cout << "z_radius " << Z_start_radius << endl;
				cout << "counter " << counter << endl;
				exit(0);
			}
		}
		double cur_distance = compute_metric_eucleidean(tunnel_map[index]->get_location_coordinates(), tunnel_map[tunnel->get_beginning_index()]->get_location_coordinates()) - tunnel_offset;
		//compute distance of current node from the tunnel beginning node
		if(DEBUG_ENABLED){
			cout << "counter " << counter << endl;
			cout << "debug" <<endl;
			cout << "tunnel_offset " << tunnel_offset << endl;
			cout << "cur distance " <<  cur_distance << endl;
			cout << "extreme point distance " << extreme_point_distance << endl;
			cout << "first node " << cutted_part_start_index << endl;
			cout << "last node" << cutted_part_end_index << endl;
			cout << "endpoit" << tunnel->get_endpoint_index() << endl;
		}
		if(TESTING_ENABLED){
			if(cur_distance > extreme_point_distance){
				cout << "TEST ERROR: node in tunnel is farther from center than extreme point in duplication check!" << endl;
				cout << "extreme point index " << extreme_point<< endl;
				cout << "endpoint_index " << cutted_part_end_index << endl;
				cout << "current index " << index << endl;
				exit(0);
			}
		}
		//assign node into interval on the basis of it's distance from beginning distance
		int interval = floor((cur_distance / extreme_point_distance) * N);
		//for the case of extreme point, which shouldn't be inserted into it's own interval(it can be any node, not just the endpoint)
		if(interval == N){
			interval--;
		}
		//std::cout << "cur " << cur_distance << " max " << extreme_point_distance << std::endl;

		//adds current nodes coordinates into vector sum representing center of gravity of current node's interval
		if(DEBUG_ENABLED){
			cout << "test2" << endl; 
			cout << "interval " << interval << endl;
			cout << "interval unfloored " << (cur_distance / extreme_point_distance) * N << endl;
			cout << "N " << N << endl;
		}
		
		if(N_counts[interval] == 0){
			delete [] N_representation[interval];
			N_representation[interval] = copy_vector(tunnel_map[index]->get_location_coordinates());
		} else {
			double* new_coordinate = add_vectors(tunnel_map[index]->get_location_coordinates(), N_representation[interval], ADDITION);
			delete [] N_representation[interval]; N_representation[interval] = new_coordinate;
		}
		N_counts[interval]++;



		
		if(index == cutted_part_end_index){
			break;
		}
		index = tunnel->get_child_index(index);
		//std::cout << "size " << tunnel_map.size() << std::endl;
		//std::cout << "N_rep loop3 " << counter << std::endl; 
		counter++; 
	}

	index = cutted_part_start_index;

	counter = 0;
	index = 0;

	//divides the sum of node position vectors in intervals by their number, so that it represents the actual center of gravity it's their interval
	while(index < N){
		if(N_counts[index] == 0){
			std::cout << "zero interval count!!" << std::endl;
			return false;
		}
		divide_vector(N_representation[index], (double) N_counts[index]);
		tunnel->add_N_point(N_representation[index]);
		//std::cout << "N_rep loop4 " << counter << std::endl; 
		counter++;
		index++;
	}

	return true;
}

//compares N point representations of two tunnels
double compute_intertunnel_distance_by_N_representations(shared_ptr<Path> path1, shared_ptr<Path> path2){

	double distance = 0;
	for(int i = 0; i < N; i++){
		//std::cout << "N " << N << std::endl; 
		distance += compute_metric_eucleidean(path1->get_N_point(i), path2->get_N_point(i));
	}
	return distance;
}


//Computes number of intervals, finds the tunnel with smallest number of nodes, that tunnel is used as basis which is divided by the INTERVAL_RATIO_CONSTANT parameter
//second argument is for situation where we don't have initialised paths yet
void decide_N(std::vector<shared_ptr<Path>>& paths, shared_ptr<Path> path){
	if(paths.size() == 0){
		N = path->get_size() / INTERVAL_RATIO_CONSTANT;
	}
	else{
		N = INT_MAX;
		for(std::vector<shared_ptr<Path>>::iterator iterator = paths.begin(); iterator != paths.end(); iterator++){ //finds path with smallest count of "atoms"
			(*iterator)->get_size() < N ? N = (*iterator)->get_size() : N = N;
		}
		N /= INTERVAL_RATIO_CONSTANT;
	}
}

//obvious
bool recompute_N_representation(std::vector<shared_ptr<Path>> existing_paths, shared_ptr<Path> path, bool try_without_changing_N){
	bool N_changed = false;
	bool successful = false;

	if(try_without_changing_N){
		path->reset_N_representation();
		successful = create_N_representation(path);
	} 
	
	while(!successful){
		N_changed = true;
		INTERVAL_RATIO_CONSTANT *= 1.1;
		decide_N(existing_paths, path);
		path->reset_N_representation();
		successful = create_N_representation(path);
	}

	return N_changed;	
}


//creates N_representation for current tunnel, it is assumed that false return value for create_N_represention means that empty interval was found, so the dividing constant is increased so that 
//we get smaller number of intervals.
void N_representation(std::vector<shared_ptr<Path>> existing_paths, shared_ptr<Path> checked_path){
	int counter = 0;
	bool was_successful = true;
	bool recompute = false;

	//find working N for the new path
	if(!create_N_representation(checked_path)){
		recompute = true;
		recompute_N_representation(existing_paths, checked_path, false);
	}

	bool recompute_checked_path = false;
	if(recompute && existing_paths.size() > 0){
		for(int i = 0; i < existing_paths.size(); i++){
			if(recompute_checked_path){
				recompute_N_representation(existing_paths, checked_path, false);
				recompute_checked_path = false;
			}
			//N changed, must start from beginning to keep sizes consistent
			if(recompute_N_representation(existing_paths, existing_paths[i], true))
			{
//				cout << "RESTART!" << endl;
				i = -1;
				recompute_checked_path = true;
				continue;
			}
//			cout << "new size " << existing_paths[i]->N_representation.size() << endl;
		}
	}
}

int is_tunnel_duplicated(shared_ptr<Path> checked_path, std::vector<shared_ptr<Path>> existing_paths){
	if(existing_paths.size() == 0){
		decide_N(existing_paths, checked_path);
		N_representation(existing_paths, checked_path);
//		cout << "N " << N << std::endl;
//		cout << "first " << checked_path->N_representation.size() << endl;
		return -1;
	}

	N_representation(existing_paths, checked_path);

//	cout << "N " << N << std::endl;
//	cout << "checked N repr " << checked_path->N_representation.size() << endl;
	for(int i = 0; i < existing_paths.size(); i++){
//		cout << "checked against N repr " << existing_paths[i]->N_representation.size() << endl;                         
		double distance = compute_intertunnel_distance_by_N_representations(existing_paths[i], checked_path);
		//std::cout << "distance " << distance << std::endl;
		//std::cout << "final N " << N << std::endl;
		if(distance < MIN_VALID_INTERTUNNEL_DISTANCE * (N / 2)){
			cout << endl;
			return i;
		}
	}
	cout << endl; 
	return -1;
}
