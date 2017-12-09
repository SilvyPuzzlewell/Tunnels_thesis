import os
import numpy as np

stats_file_path = "/home/ron/Bak/dynamic_tunnels/analysis/stats.txt"
output_file_path = "/home/ron/Bak/dynamic_tunnels/analysis/processed_stats.txt"

stats_file = open(stats_file_path, "r")
output_file = open(output_file_path, "w")

#
num_of_caver_tunnels = len(os.listdir("/home/ron/Bak/dynamic_tunnels/analysis/caver_tunnels"))
print(num_of_caver_tunnels)
iterations_number = 0
caver_tunnels_count = np.zeros(num_of_caver_tunnels)
caver_tunnels_found_in_iteration = np.zeros(num_of_caver_tunnels)
#get data
for line in stats_file:
	#1 = iteration no.? 3 = index or rrt tunnel 6 = index of related caver tunnel
	statements = line.split()
	#marks caver tunnel as found in this iteration, beginning index of caver tunnels is 1, therefore one must be subtracted to map into np array
	#if condition used to detect more rrt tunnels realating to the same caver tunnel in current iteration and therefore for 
	#disbanding possible duplicate tunnels
	if caver_tunnels_found_in_iteration[int(statements[6]) - 1] != int(statements[1]):
		caver_tunnels_count[int(statements[6]) - 1] += 1
	else:
		print("duplicated!!\n")
	caver_tunnels_found_in_iteration[int(statements[6]) - 1] = int(statements[1])
	#"last, highest will survive"
	iterations_number = int(statements[1])

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



