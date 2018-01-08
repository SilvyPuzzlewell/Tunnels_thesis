#do yourself - copy caver tunnels in analysis, choose the rest of the config file



FILE=config.txt

for number in 0.9 1.0 1.1
do
cd ~/bak_repository/dynamic_tunnels
rm ~/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '3s/.*/probe_radius '$number'/' $FILE
./run_rrt.sh 

cd ~/bak_repository/analysis
./copy_current_tunnels.sh
./analyse.sh ITERATIONS_TEST/1TQN/$number
done






