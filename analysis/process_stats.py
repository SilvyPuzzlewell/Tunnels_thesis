import os
import numpy as np
import sys
import datetime

base_path = "/home/fif/bak_repository/analysis/"
stats_file_path = base_path + "stats.txt"
pre_time_log = base_path + "found_tunnels/date_pre.log"
post_time_log = base_path + "found_tunnels/date_post.log"
time_log = base_path + "found_tunnels/time_log.log"

output_file_path = base_path + "processed_stats.txt"
output_time_file_path = base_path + "runtime.log"
output_tunnel_count_individual_iteration_file = base_path + "iterations_tunnel_counts.txt"
output_individual_iteration_time_file_path = base_path + "times.log"
if(len(sys.argv) > 1):
	output_file_path = base_path + sys.argv[1]
caver_tunnels_path = base_path + "caver_tunnels"

stats_file = open(stats_file_path, "r")
output_file = open(output_file_path, "w")
output_time_file = open(output_time_file_path, "w")
output_tunnel_count_individual_iteration_file = open(output_tunnel_count_individual_iteration_file, "w")

#
num_of_caver_tunnels = len(os.listdir(caver_tunnels_path))
print(num_of_caver_tunnels)
iterations_number = 0
caver_tunnels_count = np.zeros(num_of_caver_tunnels)
caver_tunnels_found_in_iteration = np.zeros(num_of_caver_tunnels)
caver_tunnels_found_in_iteration -= 1
#print(caver_tunnels_found_in_iteration)
#get data
iteration_tunnel_counter = 0
previous_iteration_check = 0
for line in stats_file:
	#1 = iteration no.? 3 = index or rrt tunnel 6 = index of related caver tunnel
	statements = line.split()
	iterations_number = int(statements[1])
	caver_tunnel_number = int(statements[6])
	#marks caver tunnel as found in this iteration, beginning index of caver tunnels is 1, therefore one must be subtracted to map into np array
	#if condition used to detect more rrt tunnels realating to the same caver tunnel in current iteration and therefore for 
	#disbanding possible duplicate tunnels
	if caver_tunnels_found_in_iteration[caver_tunnel_number - 1] != iterations_number:
		caver_tunnels_count[caver_tunnel_number - 1] += 1
		print(str(iterations_number)+" "+str(caver_tunnel_number))
		iteration_tunnel_counter += 1
		if iterations_number != previous_iteration_check:
			output_tunnel_count_individual_iteration_file.write(str(previous_iteration_check) + " " + str(iteration_tunnel_counter) + " " + str((iteration_tunnel_counter / num_of_caver_tunnels) * 100) + "%\n")
			iteration_tunnel_counter = 1 
		previous_iteration_check = iterations_number 
	#else:
		#print("duplicated!!\n")
	caver_tunnels_found_in_iteration[caver_tunnel_number - 1] = iterations_number
	#"last, highest will survive"
output_tunnel_count_individual_iteration_file.write(str(previous_iteration_check) + " " + str(iteration_tunnel_counter) + " " + str((iteration_tunnel_counter / num_of_caver_tunnels) * 100) + "%\n")
#process data
iterations_count = iterations_number + 1
percentage = caver_tunnels_count / iterations_count

found_tunnels_sum = 0
for n in range(0, len(caver_tunnels_count)):
	found_tunnels_sum += caver_tunnels_count[n]
	output_file.write("tunnel_num " + str(n + 1) + " is found in " + format(percentage[n] * 100, '.2f') + " percent iterations" + "\n")

found_tunnels_overall_average = found_tunnels_sum / (iterations_count * num_of_caver_tunnels)
output_file.write("overall success of rrt algoritm in finding tunnels found by caver: " + str(found_tunnels_overall_average))
output_file.close()

def compute_time_difference(time_from, time_to):

	strings_from = time_from.split()
	strings_to = time_to.split()

	dates_from = strings_from[3].split(':')
	dates_to = strings_to[3].split(':')

	hours = int(dates_to[0]) - int(dates_from[0])
	minutes = int(dates_to[1]) - int(dates_from[1])
	seconds = int(dates_to[2]) - int(dates_from[2])
	return hours*3600 + minutes*60 + seconds



pre_log = open(pre_time_log, "r")
post_log = open(pre_time_log, "r")

pre_log_line = pre_log.readlines()
pre_log_strings = pre_log_line[0].split()
pre_log_dates = pre_log_strings[3].split(':')

post_log = open(post_time_log, "r")
post_log_line = post_log.readlines()
post_log_strings = post_log_line[0].split()
post_log_dates = post_log_strings[3].split(':')

hours = int(post_log_dates[0]) - int(pre_log_dates[0])
minutes = int(post_log_dates[1]) - int(pre_log_dates[1])
seconds = int(post_log_dates[2]) - int(pre_log_dates[2])

output_time_file.write(str(hours*3600 + minutes*60 + seconds))
print(str(hours*3600 + minutes*60 + seconds))

time_log_file = open(time_log, "r")
output_individual_iteration_time_file = open(output_individual_iteration_time_file_path, "w")
time_log_lines = time_log_file.readlines()
counter = 0
while counter < len(time_log_lines):
	output_individual_iteration_time_file.write(str(compute_time_difference(time_log_lines[counter], time_log_lines[counter + 1])) + "\n")
	counter += 2










