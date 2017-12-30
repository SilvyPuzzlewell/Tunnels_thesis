#ifndef FILES
#define FILES

#include <memory>
#include <string>
#include <iostream>
#include <fstream>

#include "world.h" 

void write_paths_to_pdbs(std::string filename, std::vector<shared_ptr<Path>>& paths);
void write_to_pdb(ofstream& file, std::map<int, shared_ptr<vertex>>& atoms);
void write_clusters_to_pdbs(std::string filename, std::vector<std::vector<shared_ptr<Path>>>& clusters);

#endif