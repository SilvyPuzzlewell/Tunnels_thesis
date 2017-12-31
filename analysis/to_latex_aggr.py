import numpy
import math
from tabulate import tabulate
import sys

files_directory = "/home/fif/bak_repository/analysis/results/static/"
cur_directory = "SAMPLING_BIAS_TEST/1TQN09/"
reused = "REUSED/"
reseted = "RESETED/"
reseted_directory = files_directory + cur_directory + reseted
reused_directory = files_directory + cur_directory + reused
planFile = "processed_stats_norm.txt"



reused00_table = open(reused_directory + "00/" + planFile, "r")
reused00_lines = reused00_table.readlines()

reseted00_table = open(reseted_directory + "00/" + planFile, "r")
reseted00_lines = reseted00_table.readlines()

reused25_table = open(reused_directory + "25/" + planFile, "r")
reused25_lines = reused25_table.readlines()
#print(reused25_table.readline())
#print(reused25_table.readline())

reseted25_table = open(reseted_directory + "25/" + planFile, "r")
reseted25_lines = reseted25_table.readlines()

reused50_table = open(reused_directory + "50/" + planFile, "r")
reused50_lines = reused50_table.readlines()

reseted50_table = open(reseted_directory + "50/" + planFile, "r")
reseted50_lines = reseted50_table.readlines()

reused75_table = open(reused_directory + "75/" + planFile, "r")
reused75_lines = reused75_table.readlines()

reseted75_table = open(reseted_directory + "75/" + planFile, "r")
reseted75_lines = reseted75_table.readlines()

reused90_table = open(reused_directory + "90/" + planFile, "r")
reused90_lines = reused90_table.readlines()

reseted90_table = open(reseted_directory + "90/" + planFile, "r")
reseted90_lines = reseted90_table.readlines()

table_header = [' caver related tunnel number ', 'reused ', 'reseted', 'reused ', 'reseted', 'reused ', 'reseted', 'reused ', 'reseted', 'reused ', 'reseted']

data_table = []
data_table.append(table_header)

#-----
caver_tunnels_directory = "/home/fif/bak_repository/analysis/caver_tunnels"
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

#----



def convert_line(line):
	percentage = (float(line[5]) * 100)
	return format(percentage, '.2f') + "%"

def convert_last_line(line):
	percentage = (float(line[11]) * 100)
	return format(percentage, '.2f') + "%"

counter = 0
for line in reused00_lines:
	line = line.split()
	print(counter)
	if line[0] == "tunnel_num":
		if float(line[5]) != 0:
			caver_tunnel_length = compute_tunnel_length(line[1])
			data_line = [line[1] + " ( " + format(caver_tunnel_length, '.2f') + " ) ", 
			convert_line(reused00_lines[counter].split()), 
			convert_line(reseted00_lines[counter].split()), 
			convert_line(reused25_lines[counter].split()),
			convert_line(reseted25_lines[counter].split()),
			convert_line(reused50_lines[counter].split()), 
			convert_line(reseted50_lines[counter].split()), 
			convert_line(reused75_lines[counter].split()),
			convert_line(reseted75_lines[counter].split()),
			convert_line(reused90_lines[counter].split()),
			convert_line(reseted90_lines[counter].split())
			]
			data_table.append(data_line)
	else:
		percentage = (float(line[11]) * 100)
		data_line = [line[0] + " " + line[1],
		convert_last_line(reused00_lines[counter].split()), 
		convert_last_line(reseted00_lines[counter].split()), 
		convert_last_line(reused25_lines[counter].split()),
		convert_last_line(reseted25_lines[counter].split()),
		convert_last_line(reused50_lines[counter].split()), 
		convert_last_line(reseted50_lines[counter].split()), 
		convert_last_line(reused75_lines[counter].split()),
		convert_last_line(reseted75_lines[counter].split()),
		convert_last_line(reused90_lines[counter].split()),
		convert_last_line(reseted90_lines[counter].split()),
		]
		data_table.append(data_line)
	counter += 1
	

print(tabulate(data_table, tablefmt="latex", floatfmt = ".2f"))