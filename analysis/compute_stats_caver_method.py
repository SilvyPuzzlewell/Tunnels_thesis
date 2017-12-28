import numpy
import math 
import os
import sys
from glob import glob

#---parameters
#path to the program used in accessing necessary files
main_program_path = "/home/fif/bak_repository/analysis/"

#here must be the tunnels found by caver be copied
caver_tunnels_path = main_program_path + "caver_tunnels"

#here must be the directories in which to look for tunnels to compare to the caver ones
found_tunnels_path = main_program_path + "found_tunnels/*/"

N = 50 #default number of intervals
interval_ratio_constant = 10 #into how much intervals are the paths divided, based on maximum number of tunnel atoms in current run, example for shortest path with 1000 atoms, run will be divided into 100 intervalsw with value 10 
Z_start_radius = 2
Z_end_radius = 2
#---parameters

caver_tunnels = os.listdir(caver_tunnels_path)


def compute_distance(array1, array2):
	d1 = array1 - array2;
	d2 = d1 * d1
	d3 = d2[0] + d2[1] + d2[2]
	distance = math.sqrt(d3)
	return distance

#line is line of pdb file
def line_to_numpy(line):
	my_data = line.split()
	my_data = my_data[6:9]
	my_data = numpy.asarray(my_data, dtype='f8')
	return my_data



#tunnel is pdb file splitted to lines by readlines()
def tunnel_to_coordinate_representation(tunnel):
	tunnel_node_coordinates = []
	for tunnel_line in tunnel:
		tunnel_node_coordinates.append(line_to_numpy(tunnel_line))
	return tunnel_node_coordinates

def compute_length(tunnel):
	prev_coordinate = tunnel[0]
	length = 0
	for i in range (1, len(tunnel)):
		cur_coordinate = tunnel[i]
		length = length + compute_distance(cur_coordinate, prev_coordinate)
		prev_coordinate = cur_coordinate
	return length
def check_if_in_sphere(coordinate ,sphere_radius, sphere_coordinate):
	distance = compute_distance(coordinate, sphere_coordinate)
	if distance < sphere_radius:
		return True
	else:
		return False
def find_extreme_point(tunnel, center):
	max_distance = -1
	for coordinate in tunnel:
		cur_distance = compute_distance(coordinate, center)
		if cur_distance > max_distance:
			max_distance = cur_distance
			extreme_point = coordinate
	if max_distance <= 0:
		print("extreme point search failed!")
		sys.exit(0)

	return extreme_point

#tunnel - numpy representation
def get_N_representation(tunnel, center):
    #delete z_end interval
	index = len(tunnel) - 1
	tunnel_end_coordinate = tunnel[len(tunnel) - 1]
	tunnel_beginning_coordinate = tunnel[0]
	for i in range(0, len(tunnel)):
		if check_if_in_sphere(tunnel[index], Z_end_radius, tunnel_end_coordinate):
			tunnel = numpy.delete(tunnel, index, 0)
			index = index - 1
		else:
			break
	index = 0
	#delete z_start interval
	for coordinate in tunnel:
		if check_if_in_sphere(coordinate, Z_start_radius , tunnel_beginning_coordinate):
			tunnel = numpy.delete(tunnel, index, 0) #index is not changed because we deleting from first position
		else:
			break
	beginning_point = tunnel[0]
	extreme_point = find_extreme_point(tunnel, beginning_point)
	extreme_point_distance = compute_distance(extreme_point, beginning_point)

	intervals = numpy.zeros((N, 3))
	counts = numpy.zeros(N)

	#assign points to intervals  
	for coordinate in tunnel:
		distance = compute_distance(coordinate, beginning_point)
		#print("------")
		interval = int(numpy.floor((distance / extreme_point_distance) * N))
		#print(interval)
		#print(len(counts))
		#extreme point element would belong to new interval..
		if interval == len(counts):
			interval = interval - 1

		counts[interval] = counts[interval] + 1
		intervals[interval] = intervals[interval] + coordinate
	for i in range(0, len(intervals)):
		if counts[i] == 0:
			print("zero interval count!!")
			exit(0)
		intervals[i] = intervals[i] / counts[i]
	return intervals

def min_number_of_atoms(directory_path, directory_content):
	min_count = sys.maxsize
	for tunnel in directory_content:
		num_lines = sum(1 for line in open(directory_path + "/" + tunnel))
		#print(num_lines)
		if num_lines < min_count:
			min_count = num_lines
	return min_count

def inter_tunnel_representation_distance(tunnel_N1, tunnel_N2):
	distance = 0
	for n in range(0, N):
		distance = distance + ((n + 1) / N) * compute_distance(tunnel_N1[n], tunnel_N2[n])
	return distance



#x = numpy.array([2, 6, 8])
#y = numpy.array([3, 6, 8])
#print(x - y)

#used only as indicator of number of directories to process, loads every directory into list
def load_directories_list():
	return (glob(found_tunnels_path))

output_file = "stats.txt"
iteration = 0
output = open(output_file, "w")

#at this time, only the beginning vertex is used as the center of gravity of the whole thing
config_file = open(main_program_path + "config.txt")
config = config_file.readlines()
center_coordinates = numpy.asarray([float((config[4].split()[1])), float((config[5].split()[1])), float((config[6].split()[1]))])
config_file.close()

coords = tunnel_to_coordinate_representation((open(caver_tunnels_path + "/1").readlines()))

def compute_tunnel_length_width(tunnel):
	counter = 0
	length = 0
	for coordinate in tunnel:
		if counter != 0:
			length = length + compute_distance(tunnel[counter], tunnel[counter - 1])
		counter = counter + 1
	return length


caver_tunnels = os.listdir(caver_tunnels_path)
N = int(numpy.floor(min_number_of_atoms(caver_tunnels_path, caver_tunnels) / interval_ratio_constant))
print(N)
caver_tunnel_N_representations_list = []
indices_mapping = []
lengths = open("lengths.txt", "w")
for tunnel in caver_tunnels:
	caver_tunnel = open(caver_tunnels_path + "/" + tunnel, "r")
	#print(tunnel)
	tunnel_array = tunnel_to_coordinate_representation(caver_tunnel.readlines())
	tunnel_representation = get_N_representation(tunnel_array, center_coordinates)
	caver_tunnel_N_representations_list.append(tunnel_representation)
	indices_mapping.append(int(tunnel))
	lengths.write("tunnel " + tunnel + ", length: " + str(compute_length(tunnel_array)) + "\n")
	#print(indices_mapping)




iterations = (glob(found_tunnels_path))
output_file = "stats.txt"
iteration = 0
output = open(output_file, "w")

for iteration_c in iterations: #iterating through folders containing tunnels found by individual runs of tunnel detection algorithm
	counter = 0

	my_tunnels_path = main_program_path + "found_tunnels/tunnels" + str(iteration)
	my_tunnels = os.listdir(my_tunnels_path)

	for file in my_tunnels: #iterating through individual tunnels found in one run of tunnel detection algorithm
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

		tunnel_distance = sys.maxsize

		related_tunnel_index = -1

		tunnel_array = tunnel_to_coordinate_representation(my_tunnel.readlines())
		print(file)
		tunnel_representation = get_N_representation(tunnel_array, center_coordinates)

		for n in range(0, len(caver_tunnel_N_representations_list)):
			cur_distance = inter_tunnel_representation_distance(tunnel_representation, caver_tunnel_N_representations_list[n])
			#print(cur_distance)
			if cur_distance < tunnel_distance:
				tunnel_distance = cur_distance
				related_tunnel_index = indices_mapping[n]
        #outputting relation of current tunnel to most similiar tunnel found by caver
		#print(related_tunnel_index)
		#print(related_tunnel_index)
		output_line = "iteration " + str(iteration) + " my_tunnel " + str(file_number) + " related to " + str(related_tunnel_index) + " length " + str(compute_length(tunnel_array)) + " distance " + str(tunnel_distance) + "\n"
		print(output_line)
		output.write(output_line)
		counter += 1 #counter of tunnels found by rrt in current run
	iteration += 1 #counter of rrt runs to be analyzed
output.close()

