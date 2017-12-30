cd clusters/clusters0
directories=$(find * -type d)

for directory in $directories;
do
	echo $directory
	cd $directory
	word_count=$(ls | wc -w)
	for i in `seq 2 $word_count`;
	do
		cat "trajectory$i.pdb" >> trajectory1.pdb
		rm "trajectory$i.pdb"
	done
	cd ..
done



