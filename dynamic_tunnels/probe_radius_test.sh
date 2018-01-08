#$1 - protein



FILE=config.txt
PATH_TO_CAVER_TUNNELS=caver_tunnels_database/$1/

for number in 0.6 0.7 0.8 0.9 1.0 1.1
do
cd ~/bak_repository/dynamic_tunnels
rm ~/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '3s/.*/probe_radius '$number'/' $FILE
./run_rrt.sh 

cd ~/bak_repository/analysis
rm caver_tunnels/*
cp $PATH_TO_CAVER_TUNNELS/$number/* caver_tunnels/

./copy_current_tunnels.sh
./analyse.sh ITERATIONS_TEST/1TQN/$number
done
