import numpy
import math
from tabulate import tabulate
import sys

if len(sys.argv) > 1:
	planFile = sys.argv[1]
else:
	planFile = "processed_stats.txt"

file_table = open(planFile)
file_lines = file_table.readlines()


table_header = [' caver related tunnel number ', ' found % times ']

data_table = []
data_table.append(table_header)

for line in file_lines:
	line = line.split()
	if line[0] == "tunnel_num":
		data_line = [line[1], str(math.floor(float(line[5]) * 100)) + "%"]
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