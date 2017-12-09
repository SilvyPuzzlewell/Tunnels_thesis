import os

#edited tunnels are saved in ./caver_tunnels folder named by their respective number
newpath = '/home/ron/Bak/dynamic_tunnels/analysis/caver_tunnels' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
tunnels_directory = os.listdir(newpath)
for file in tunnels_directory:
	os.remove(newpath+"/"+file)
#path to caver tunnels
dir_path = "/home/ron/Bak/caver/examples/QUICK_START/inputs/out/data/clusters_timeless"
tunnels = os.listdir(dir_path)

for file in tunnels:
	beginning = file[:4]
	if beginning == "tun_":
		file_number = int(file[7:10])
		Out_file = "caver_tunnels/" + str(file_number)
		data = open(dir_path + "/" + file, "r")
		output = open(Out_file, "w")
		for line in data:
			words = line.split()
			if len(words) != 0 and words[0] == "ATOM":
				output.write(line)
		data.close()
		output.close()

