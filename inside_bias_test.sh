#$1 - protein



FILE=config.txt
PATH_TO_CAVER_TUNNELS=caver_tunnels_database/$1/

sed -i '13s/.*/reseted_tree_mode 1/' $FILE

for number in 0 0.25 0.5 0.75 0.9 0.99
do
cd ~/bak_repository/dynamic_tunnels_copy#1
rm ~/bak_repository/dynamic_tunnels_copy#1/runtime_stats.txt
sed -i '11s/.*/probe_radius '$number'/' $FILE
./run_rrt.sh 

cd ~/bak_repository/analysis#1
./copy_current_tunnels.sh
./analyse.sh ITERATIONS_TEST/1TQN/RESETED/$number
done
