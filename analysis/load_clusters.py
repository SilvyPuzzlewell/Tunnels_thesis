from pymol import cmd
import os 
from glob import glob

dir_path = "/home/fif/bak_repository/analysis/clusters/clusters0"
cluster_directories = glob(dir_path + "/*/")
colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
           'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', \
           'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',    \
           'wheat', 'white', 'grey' ]
counter = 0

for cluster in cluster_directories:
	tunnels = os.listdir(cluster)
	cmd.cd(cluster)
	print(cluster)
	for x in tunnels:
		select = x[:-4]
		cmd.load(x)

		print(x)
		print(select)

		
		cmd.color(colours[counter], select)
		cmd.alter(select, 'vdw=b')
		cmd.show("spheres", select)
		cmd.set_name(select, select + "#" +str(counter))
	counter = counter + 1
	if counter > len(colours) - 1:
		counter = 0

cmd.cd("/home/fif/bak_repository/analysis")