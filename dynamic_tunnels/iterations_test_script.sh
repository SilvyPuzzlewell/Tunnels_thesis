#do yourself - copy caver tunnels in analysis, choose the rest of the config file



FILE=config.txt

for number in 20000 50000 100000 200000 500000 1000000 
do
cd ~/bak_repository/dynamic_tunnels
rm ~/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '2s/.*/iterations '$number'/' $FILE
./run_rrt.sh 

cd ~/bak_repository/analysis
./copy_current_tunnels.sh
./analyse.sh ITERATIONS_TEST/$1/$number
done






