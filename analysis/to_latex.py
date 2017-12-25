import numpy
import math
from tabulate import tabulate
import sys

caver_tunnels_directory = "/home/fif/bak_repository/analysis/caver_tunnels"

if len(sys.argv) > 1:
	planFile = sys.argv[1]
else:
	planFile = "processed_stats.txt"

file_table = open(planFile)
file_lines = file_table.readlines()

table_header = [' caver related tunnel number ', ' found % times ']

data_table = []
data_table.append(table_header)

def line_to_numpy(line):
	my_data = line.split()
	my_data = my_data[6:9]
	my_data = numpy.asarray(my_data, dtype='f8')
	return my_data

def tunnel_to_coordinate_representation(tunnel):
	tunnel_node_coordinates = []
	for tunnel_line in tunnel:
		tunnel_node_coordinates.append(line_to_numpy(tunnel_line))
	return tunnel_node_coordinates

def compute_distance(array1, array2):
	d1 = array1 - array2;
	d2 = d1 * d1
	d3 = d2[0] + d2[1] + d2[2]
	distance = math.sqrt(d3)
	return distance

def compute_length(tunnel):

	tunnel_array = tunnel_to_coordinate_representation(tunnel.readlines())
	prev_coordinate = tunnel_array[0]
	length = 0
	for i in range (1, len(tunnel_array)):
		cur_coordinate = tunnel_array[i]
		length = length + compute_distance(cur_coordinate, prev_coordinate)
		prev_coordinate = cur_coordinate
	return length

def compute_tunnel_length(file):
	caver_tunnel = open(caver_tunnels_directory + "/" + file, "r")
	return compute_length(caver_tunnel)


for line in file_lines:
	line = line.split()
	if line[0] == "tunnel_num":
		if float(line[5]) != 0:
			caver_tunnel_length = compute_tunnel_length(line[1])
			data_line = [line[1] + " ( " + format(caver_tunnel_length, '.2f') + " ) ", str(math.floor(float(line[5]) * 100)) + "%"]
			data_table.append(data_line)
	else:
		data_line = [line[0] + " " + line[1], str(math.floor(float(line[11]) * 100)) + "%"]
		data_table.append(data_line)
	

print(tabulate(data_table, tablefmt="latex", floatfmt = ".2f"))
print("\multicolumn{2}{c}{Item} \\")



"""for line in file_table:

table = numpy.loadtxt(planFile)
print(planFile)
"""
"""
data_normal = numpy.loadtxt(planFile_normal)
data_path = numpy.loadtxt(planFile_path)

table_header = (sys.argv[1], ' RRT ', ' RRT_Path ')

mean_normal = numpy.mean(data_normal, axis=0)
stddev_normal = numpy.std(data_normal, axis=0)

mean_path = numpy.mean(data_path, axis=0)
stddev_path = numpy.std(data_path, axis=0)

table_mat = [table_header, [' number of iterations ', "{0:.2f}".format(mean_normal[0]) + " " + "{0:.2f}".format(stddev_normal[0]), "{0:.2f}".format(mean_path[0]) + " " + "{0:.2f}".format(stddev_path[0])], \
[' time ', "{0:.2f}".format(mean_normal[1]) + " " + "{0:.2f}".format(stddev_normal[1]), "{0:.2f}".format(mean_path[1]) + " " + "{0:.2f}".format(stddev_path[1])], \
[' tree size ', "{0:.2f}".format(mean_normal[2]) + " " + "{0:.2f}".format(stddev_normal[2]), "{0:.2f}".format(mean_path[2]) + " " + "{0:.2f}".format(stddev_path[2])], \
[' distance to goal ', "{0:.2f}".format(mean_normal[3]) + " " + "{0:.2f}".format(stddev_normal[3]), "{0:.2f}".format(mean_path[3]) + " " + "{0:.2f}".format(stddev_path[3])]]
print(tabulate(table_mat, tablefmt="latex", floatfmt=".2f"))
"""
