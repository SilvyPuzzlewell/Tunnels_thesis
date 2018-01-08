import numpy
import math
from tabulate import tabulate
import sys

molecule="1MXTa"
molecule_directory = "1MXTa/"

files_directory = "/home/fif/bak_repository/analysis/results/static/"
cur_directory = "PROBE_TEST/" + molecule + "/"
caver_tunnels_directory = "/home/fif/bak_repository/analysis/caver_tunnels_database/" + molecule + "/"
reused = ""
reseted = ""
tunnel_directory = files_directory + cur_directory
caver_directory = caver_tunnels_directory + molecule_directory
planFile = "processed_stats_norm.txt"
timeFile = "runtime.log"



radius1_directory = "0.6/"
radius2_directory = "0.7/" #sudy patri caveru
radius3_directory = "0.8/"

probe1_file = open(tunnel_directory + radius1_directory + planFile, "r")
probe1_time_file = open(tunnel_directory + radius1_directory + timeFile, "r")
probe2_file = open(tunnel_directory + radius2_directory + planFile, "r")
probe2_time_file = open(tunnel_directory + radius2_directory + timeFile, "r")
probe3_file = open(tunnel_directory + radius3_directory + planFile, "r")
probe3_time_file = open(tunnel_directory + radius3_directory + timeFile, "r")



custom_style = True
another_title = True

non_caver_lines = 4
title_index = 2
style = "\\begin{tabular}{|ll|ll|ll|}"
title = "\\multicolumn{6}{|c|}{TITLE} \\\\"
Another_title = "\\multicolumn{2}{|c|}{0.6} & \\multicolumn{2}{c|}{0.7} & \\multicolumn{2}{|c|}{0.8}\\\\"
table_header = [' CAVER tunnel n. ', ' 0.6 ', ' 0.7', '0.8', '0.9', '1.0 ', '1.1']


data_table = []
data_table.append(table_header)

#-----

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
	caver_tunnel = open(caver_tunnels_directory  + file, "r")
	return compute_length(caver_tunnel)

#----



def convert_line(line_in):
	line = line_in.split()
	percentage = (float(line[5]))
	return format(percentage, '.0f') + "%"

def convert_last_line(line_in):
	line = line_in.split()
	percentage = (float(line[11])* 100)
	return format(percentage, '.2f') + "%"

def return_percentage(line, len, num_line):

	if len <= num_line:
		return 0
	line = line.split()
	return float(line[5]) 

def print_lines(lines):
	for line in lines:
		print(line)


column_tables = {probe1_file, probe2_file, probe3_file}


def return_caver_line(probe_radius_directory, file_number):

	length = compute_tunnel_length(file_number)
	return file_number[4:] + " ( " + format(length, '.2f') + " ) "

def return_check(probe_table, counter, mode):
	#print(probe_table)
	
	
	if(len(probe_table) > counter):

		return probe_table[counter][mode]
	else:
		return "" 


counter = 0
caver_lines = 0


data_table = []
#subtables

probe1_table = []
probe2_table = []
probe3_table = []
probe1_lines = probe1_file.readlines()
probe2_lines = probe2_file.readlines()
probe3_lines = probe3_file.readlines()

counter1 = 0
final = 0
for line in probe1_lines:
	if(line[:10] == "tunnel_num"):
		if(return_percentage(line ,len(probe1_lines), counter1) != 0):
			data_line = [return_caver_line(radius1_directory, "0.6/" + str(counter1 + 1)), convert_line(line)]
			probe1_table.append(data_line)
			counter1 += 1
	else:
		data_line = ["", convert_last_line(line)]
		if counter1 > final:
			final = counter1 

counter2 = 0
for line in probe2_lines :
	if(line[:10] == "tunnel_num"):
		if(return_percentage(line, len(probe2_lines), counter2) != 0):
			data_line = [return_caver_line(radius2_directory, "0.7/" + str(counter2 + 1)), convert_line(line)]
			probe2_table.append(data_line)
			counter2 += 1
	else:
		data_line = ["", convert_last_line(line)]
		if counter2 > final:
			final = counter2 



counter3 = 0
for line in probe3_lines:
	if(line[:10] == "tunnel_num"):
		if(return_percentage(line, len(probe3_lines), counter3) != 0):
			data_line = [return_caver_line(radius3_directory, "0.8/" + str(counter3 + 1)), convert_line(line)]
			probe3_table.append(data_line)
			counter3 += 1
	else:
		data_line = ["", convert_last_line(line)]
		if counter3 > final:
			final = counter3 




counter = 0
while counter < final:
	data_line = [return_check(probe1_table, counter, 0), return_check(probe1_table, counter, 1),
			return_check(probe2_table, counter, 0), return_check(probe2_table, counter, 1),
			return_check(probe3_table, counter, 0), return_check(probe3_table, counter, 1)		
			]
	counter += 1
	data_table.append(data_line)



table = tabulate(data_table, tablefmt="latex", floatfmt = ".2f")
table_lines = table.split("\n")
if custom_style:
	table_lines[0] = style
#---insert title
table_lines.insert(title_index, title)
table_lines.insert(title_index + 1, "	\hline")
non_caver_lines += 2

if another_title:
	table_lines.insert(title_index + 2, Another_title)
	table_lines.insert(title_index + 3, "	\hline")
	non_caver_lines += 2

#---

#--- insert runtime stats
after_comparison_index = non_caver_lines + final - 2

data_line = "\\multicolumn{2}{|c|}{" + str(probe1_time_file.read()) +"s} & \\multicolumn{2}{c|}{" + str(probe2_time_file.read()) +"s} & \\multicolumn{2}{|c|}{" + str(probe3_time_file.read()) +"s}\\\\"
table_lines.insert(after_comparison_index, data_line)
#----

print_lines(table_lines)


#print(tabulate(data_table, tablefmt="latex", floatfmt = ".2f"))
