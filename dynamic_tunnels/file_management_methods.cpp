#include <memory>
#include <iostream>   
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "world.h"
#include "rrt_methods.h"
#include "support_methods.h"


void write_to_pdb(ofstream& file, std::map<int, shared_ptr<vertex>>& atoms){
  std::map<int, shared_ptr<vertex>>::iterator iterator;
  for(iterator = atoms.begin(); iterator != atoms.end(); iterator++){
      file << setprecision(3) << fixed;
      file << "ATOM";

      file.width(7); file << iterator->second->get_first_frame();
      float c1 = 1;
      file << "  N   ILE E  16";
      file.width(12); file <<  iterator->second->get_location_coordinates()[0]; file.width(8); file <<  iterator->second->get_location_coordinates()[1]; file.width(8); file <<  iterator->second->get_location_coordinates()[2]; file.width(2); file.precision(2); file << "  "; file.width(4); file << c1;
      file.width(1); file << " "; file.width(5); file << iterator->second->get_radius(); file.width(12); file << "N"<<endl;      
  }
}


void write_paths_to_pdbs(std::string filename, std::vector<shared_ptr<Path>>& paths){

  int counter = 0;
  std::vector<shared_ptr<Path>>::iterator iterator;
  int i = 0;
  for(iterator = paths.begin(); iterator != paths.end(); iterator++){
    if(!test_path_noncolliding_static(*iterator)){
      cout << "path colliding before write!" <<endl;
      create_segfault();
   }
    std::ofstream file;
    std::string filename_appended(filename + to_string(i + 1));
    filename_appended.append(".pdb");
    file.open(filename_appended); 
    std::cout << i << " " << filename_appended << endl;
    write_to_pdb(file, (*iterator)->get_vertices());
    file.close();
    i++;
 }
}

void write_clusters_to_pdbs(std::string directory_template, std::vector<std::vector<shared_ptr<Path>>>& clusters){
  for(int i = 0; i < clusters.size(); i++){
    std::string directory(directory_template + to_string(i + 1));
    std::string filename(directory + "/trajectory");
    std::string system_call("mkdir " + directory);
    system(system_call.c_str());
    write_paths_to_pdbs(filename, clusters[i]);
  }
}
