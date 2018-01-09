import numpy
import math
from tabulate import tabulate
import sys

files_directory = "/home/fif/bak_repository/analysis/results/static/"
cur_directory = "INSIDE_BIAS_TEST/1MXTa/"
reused = "REUSED/"
reseted = "RESETED/"
reseted_directory = files_directory + cur_directory + reseted
reused_directory = files_directory + cur_directory + reused
planFile = "processed_stats_norm.txt"
timeFile = "runtime.log"



column1_directory = "0.75/"
column2_directory = "0.75/"
column3_directory = "0.9/"
column4_directory = "0.9/"
column5_directory = "300k/"
column6_directory = "300k/"



column1_table = open(reused_directory + column1_directory + planFile, "r")
column1_time = open(reused_directory + column1_directory + timeFile, "r")
column1_lines = column1_table.readlines()

column2_table = open(reseted_directory + column2_directory + planFile, "r")
column2_time = open(reseted_directory + column2_directory + timeFile, "r")
column2_lines = column2_table.readlines()

column3_table = open(reused_directory + column3_directory + planFile, "r")
column3_time = open(reused_directory + column3_directory + timeFile, "r")
column3_lines = column3_table.readlines()
#print(reused25_table.readline())
#print(reused25_table.readline())

column4_table = open(reseted_directory + column4_directory + planFile, "r")
column4_time = open(reseted_directory + column4_directory + timeFile, "r")
column4_lines = column4_table.readlines()

column5_table = open(reused_directory + column5_directory + planFile, "r")
column5_time = open(reused_directory + column5_directory + timeFile, "r")
column5_lines = column5_table.readlines()

column6_table = open(reseted_directory + column6_directory + planFile, "r")
column6_time = open(reseted_directory + column6_directory + timeFile, "r")
column6_lines = column6_table.readlines()




custom_style = True
another_title = True

non_caver_lines = 4
title_index = 2
style = "\\begin{tabular}{|l|ll|ll|ll|}"
title = "\\multicolumn{7}{|c|}{TITLE} \\\\"
Another_title = "\\multicolumn{1}{|c|}{} & \\multicolumn{2}{c|}{0.75} & \\multicolumn{2}{|c|}{0.9} & \\multicolumn{2}{c|}{300k}\\\\"
table_header = [' caver tunnel n. ', 'reused ', 'reseted', 'reused ', 'reseted', 'reused ', 'reseted']


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



def convert_line(line_in):
	line = line_in.split()
	percentage = (float(line[5]))
	return format(percentage, '.0f') + "%"

def convert_last_line(line_in):
	line = line_in.split()
	percentage = 0
	if len(line) > 8:
		percentage = (float(line[11])* 100)
	return format(percentage, '.2f') + "%"

def return_percentage(line_in):
	line = line_in.split()
	return float(line[5]) 

def print_lines(lines):
	for line in lines:
		print(line)


counter = 0
caver_lines = 0
for line in column1_lines:
	line = line.split()
	print(counter)
	if line[0] == "tunnel_num":

		if (return_percentage(column1_lines[counter]) != 0 or return_percentage(column2_lines[counter]) != 0 or return_percentage(column3_lines[counter]) != 0 
		or return_percentage(column4_lines[counter]) != 0 or return_percentage(column5_lines[counter]) != 0 or return_percentage(column6_lines[counter]) != 0):
			caver_tunnel_length = compute_tunnel_length(line[1])
			data_line = [line[1] + " ( " + format(caver_tunnel_length, '.2f') + " ) ", 
			convert_line(column1_lines[counter]), 
			convert_line(column2_lines[counter]), 
			convert_line(column3_lines[counter]),
			convert_line(column4_lines[counter]),
			convert_line(column5_lines[counter]), 
			convert_line(column6_lines[counter])
			]
			data_table.append(data_line)
			caver_lines += 1
	else:
		percentage = (float(line[11]) * 100)
		data_line = [line[0] + " " + line[1],
		convert_last_line(column1_lines[counter]), 
		convert_last_line(column2_lines[counter]), 
		convert_last_line(column3_lines[counter]),
		convert_last_line(column4_lines[counter]),
		convert_last_line(column5_lines[counter]), 
		convert_last_line(column6_lines[counter])
		]
		data_table.append(data_line)
	counter += 1

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
after_comparison_index = non_caver_lines + caver_lines

data_line = (" runtime " + " & " + str(column1_time.read()) + "s & " + str(column2_time.read()) + "s & " + str(column3_time.read()) + "s & " + str(column4_time.read())
+ "s & " + str(column5_time.read()) + "s & " + str(column6_time.read()) + " \\\\")
table_lines.insert(after_comparison_index, data_line)
#----

print_lines(table_lines)


#print(tabulate(data_table, tablefmt="latex", floatfmt = ".2f"))
