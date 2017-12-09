from pymol import cmd
from os import listdir

dir_path = "/home/ron/Bak/dynamic_tunnels/found_tunnels/tunnels0"
frames_count = 50
increment = 255 / frames_count

for i in range(1,frames_count+1):
	cmd.set_color(str("frame" + str(i)), [255 - increment * (i - 1), 0, 0 + increment * (i - 1)])

tunnels = os.listdir(dir_path)
cmd.cd(dir_path)

counter = 0
for x in tunnels:
	cmd.load(x)
	select = x[:-4]
	cmd.color(colours[counter], select)
	cmd.alter(select, 'vdw=b')
	cmd.show("spheres", select)
	counter = counter + 1
	if counter > len(colours) - 1:
		counter = 0
cmd.cd("/home/ron/Bak/dynamic_tunnels/")
