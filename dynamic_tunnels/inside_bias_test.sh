#$1 - protein



FILE=config.txt
PATH_TO_CAVER_TUNNELS=caver_tunnels_database/$1/

sed -i '13s/.*/reseted_tree_mode 0/' $FILE

for number in 0 0.25 0.5 0.75 0.9 0.99
do
cd ~/bak_repository/dynamic_tunnels
rm ~/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '11s/.*/inside_sampling_bias '$number'/' $FILE
./run_rrt.sh 

cd ~/bak_repository/analysis
./copy_current_tunnels.sh
./analyse.sh INSIDE_BIAS_TEST/1AKD/REUSED/$number
done

for number in 0 0.25 0.5 0.75 0.9 0.99
do
cd ~/bak_repository/dynamic_tunnels
rm ~/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '2s/.*/iterations 300000/' $FILE
sed -i '11s/.*/inside_sampling_bias 0/' $FILE
./run_rrt.sh 

cd ~/bak_repository/analysis
./copy_current_tunnels.sh
./analyse.sh INSIDE_BIAS_TEST/1AKD/REUSED/CROSSTEST/300k
