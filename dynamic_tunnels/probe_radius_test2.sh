#$1 - protein



FILE=config.txt
PATH_TO_CAVER_TUNNELS=caver_tunnels_database/$1/

for number in 0 0.25 0.5 0.75 0.90 0.99
do
cd ~/bak_repository/dynamic_tunnels
rm ~/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '11s/.*/inside_sampling_bias '$number'/' $FILE
./run_rrt.sh 

cd ~/bak_repository/analysis
./analyse.sh INSIDE_BIAS_TEST/1AKD06/$number
done


