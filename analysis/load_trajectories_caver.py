from pymol import cmd
from os import listdir

dir_path = "/home/fif/bak_repository/analysis/"
tunnels = os.listdir(dir_path)
cmd.cd(dir_path)

colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
           'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', \
           'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',    \
           'wheat', 'white', 'grey' ]
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
cmd.cd("/home/fif/bak_repository/")
