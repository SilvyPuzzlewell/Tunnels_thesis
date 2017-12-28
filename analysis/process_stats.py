import os
import numpy as np
import sys
import datetime

base_path = "/home/fif/bak_repository/analysis/"
stats_file_path = base_path + "stats.txt"
pre_time_log = base_path + "found_tunnels/date_pre.log"
post_time_log = base_path + "found_tunnels/date_post.log"

output_file_path = base_path + "processed_stats.txt"
output_time_file_path = base_path + "runtime.log"
if(len(sys.argv) > 1):
	output_file_path = base_path + sys.argv[1]
caver_tunnels_path = base_path + "caver_tunnels"

stats_file = open(stats_file_path, "r")
output_file = open(output_file_path, "w")
output_time_file = open(output_time_file_path, "w")

#
num_of_caver_tunnels = len(os.listdir(caver_tunnels_path))
print(num_of_caver_tunnels)
iterations_number = 0
caver_tunnels_count = np.zeros(num_of_caver_tunnels)
caver_tunnels_found_in_iteration = np.zeros(num_of_caver_tunnels)
caver_tunnels_found_in_iteration -= 1
#print(caver_tunnels_found_in_iteration)
#get data
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
	#else:
		#print("duplicated!!\n")
	caver_tunnels_found_in_iteration[caver_tunnel_number - 1] = iterations_number
	#"last, highest will survive"

#process data
iterations_count = iterations_number + 1
percentage = caver_tunnels_count / iterations_count

found_tunnels_sum = 0
for n in range(0, len(caver_tunnels_count)):
	found_tunnels_sum += caver_tunnels_count[n]
	output_file.write("tunnel_num " + str(n + 1) + " is found in " + str(percentage[n]) + " percent iterations" + "\n")

found_tunnels_overall_average = found_tunnels_sum / (iterations_count * num_of_caver_tunnels)
output_file.write("overall success of rrt algoritm in finding tunnels found by caver: " + str(found_tunnels_overall_average))
output_file.close()




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







