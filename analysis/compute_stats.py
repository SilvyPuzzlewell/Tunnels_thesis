import numpy
import math 
import os
import re
from scipy import spatial as sp
from glob import glob

#here must tunnels found by caver be copied
caver_tunnels_path = "/home/fif/bak_repository/analysis/caver_tunnels"

rrt_tunnels_base_path = "/home/fif/bak_repository/analysis/found_tunnels/"
rrt_tunnels_path = rrt_tunnels_base_path + "*/"
#naming converntion for tunnel directories, the prefix after which the tunnel number is
rrt_tunnels_namepattern = rrt_tunnels_base_path + "tunnels"

VALID_TUNNEL_PATTERN = re.compile("(trajectory\d+\.pdb)|(\d+)")


def compute_distance(array1, array2):
	d1 = array1 - array2;
	d2 = d1 * d1
	d3 = d2[0] + d2[1] + d2[2]
	distance = math.sqrt(d3)
	return distance

def line_to_numpy(line):
	my_data = line.split()
	my_data = my_data[6:9]
	my_data = numpy.asarray(my_data, dtype='f8')
	return my_data

#beginning sequence is used to identify correct files for kd tree construction
def init_kd_trees(path, directory):
	kd_trees = []
	tree_data = []
	caver_counter = 0


	for tunnel in directory:
		data = open(path + "/" + directory[caver_counter], "r")
		if not VALID_TUNNEL_PATTERN.match(directory[caver_counter]):
			data.close()
			caver_counter += 1
			continue
		for tunnel_line in data:						
			tunnel_data = line_to_numpy(tunnel_line)
			tree_data.append(tunnel_data)
		search_tree = sp.cKDTree(tree_data, copy_data="true")
		kd_trees.append(search_tree)
		data.close()
		tree_data = []
		caver_counter += 1
	return  kd_trees

#computes one way distance using kd tree of second input file to find second file node with shortest distance belonging to first file, out[0] = distance between "files", out[1] = length of tunnel in file
def file_to_file_distance(tunnel_file, kd_tree):
	counter = 0
	tunnel_distance = 0
	tunnel_length = 0
	for tunnel_line in tunnel_file:
		data = line_to_numpy(tunnel_line)
		if counter > 0:
			tunnel_distance += kd_tree.query(data)[0]
			tunnel_length += compute_distance(prev_data, data)
		prev_data = data
		counter = counter + 1
	ret = []
	ret.append(tunnel_distance)
	ret.append(tunnel_length)
	return ret




caver_tunnels = os.listdir(caver_tunnels_path)


#for nearest neighbor search when checking nearest neighbor of my tunnels in caver tunnel, caver tunnels are still the same, therefore no need for reconstructing
caver_kd_trees = init_kd_trees(caver_tunnels_path, caver_tunnels)

#used only as indicator of number of directories to process, loads every directory into list
iterations = (glob(rrt_tunnels_path))

output_file = "stats.txt"
iteration = 0
output = open(output_file, "w")
for iteration_c in iterations: #iterating through folders containing tunnels found by individual runs of tunnel detection algorithm
	counter = 0

	my_tunnels_path = rrt_tunnels_namepattern + str(iteration)
	my_tunnels = os.listdir(my_tunnels_path)
	my_kd_trees = init_kd_trees(my_tunnels_path, my_tunnels)

	for file in my_tunnels: #iterating through individual tunnels found in one run of tunnel detection algorithm
		if not VALID_TUNNEL_PATTERN.match(file):
			print("invalid file " + file)
			continue
		beginning = file[:10]

		file_number = ""
		if beginning == "trajectory": #to print in stats
			file_end = file[10:]
			for c in file_end:
				if c == ".":
					break
				file_number = file_number + c
		else:
			continue
		
		my_tunnel = open(my_tunnels_path + "/" + file, "r")

		tunnel_distance = 0
		tunnel_length = 0
		ratio = float("inf")
		related_tunnel_index = -1

		print(file)
		for n in range(0, len(caver_kd_trees)): #finding shortest distance in of current tunnel by iterating through individual tunnels found by caver algorithm
			dist = file_to_file_distance(my_tunnel, caver_kd_trees[n])
			caver_tunnel = open(caver_tunnels_path + "/" + caver_tunnels[n], "r")
			print("here " + str(counter))
			print(len(my_kd_trees))
			dist2 = file_to_file_distance(caver_tunnel, my_kd_trees[counter])
			caver_tunnel.close()
			
			my_tunnel.seek(0)

			cur_ratio = dist[0] / dist[1]
			cur_ratio2 = dist2[0] / dist2[1]
			if cur_ratio2 > cur_ratio:
				cur_ratio = cur_ratio2

			if cur_ratio < ratio:
				ratio = cur_ratio
				related_tunnel_index = n
			tunnel_length_out = dist[1]
        #outputting relation of current tunnel to most similiar tunnel found by caver
		print(related_tunnel_index)
		output_line = "iteration " + str(iteration) + " my_tunnel " + str(file_number) + " related to " + caver_tunnels[related_tunnel_index] + " length " + str(tunnel_length_out) + " ratio " + str(ratio) + "\n"
		print(output_line)
		output.write(output_line)
		counter += 1 #counter of tunnels found by rrt in current run
	iteration += 1 #counter of rrt runs to be analyzed
output.close()







		
		
