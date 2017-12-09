



FILE=/home/ron/Bak/dynamic_tunnels/config.txt
rm /home/ron/Bak/dynamic_tunnels/runtime_stats.txt
sed -i '2s/.*/iterations 20000/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_20kiterations_3SB
mv found_tunnels/* computed_at_night/1TQN_09_20kiterations_3SB
mv runtime_stats.txt computed_at_night/1TQN_09_20kiterations_3SB


FILE=/home/ron/Bak/dynamic_tunnels/config.txt
rm /home/ron/Bak/dynamic_tunnels/runtime_stats.txt
sed -i '2s/.*/iterations 50000/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_50kiterations_3SB
mv found_tunnels/* computed_at_night/1TQN_09_50kiterations_3SB
mv runtime_stats.txt computed_at_night/1TQN_09_50kiterations_3SB


FILE=/home/ron/Bak/dynamic_tunnels/config.txt
rm /home/ron/Bak/dynamic_tunnels/runtime_stats.txt
sed -i '2s/.*/iterations 100000/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_100kiterations_3SB
mv found_tunnels/* computed_at_night/1TQN_09_100kiterations_3SB
mv runtime_stats.txt computed_at_night/1TQN_09_100kiterations_3SB

FILE=/home/ron/Bak/dynamic_tunnels/config.txt
rm /home/ron/Bak/dynamic_tunnels/runtime_stats.txt
sed -i '2s/.*/iterations 200000/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_200kiterations_3SB
mv found_tunnels/* computed_at_night/1TQN_09_200kiterations_3SB
mv runtime_stats.txt computed_at_night/1TQN_09_200kiterations_3SB


FILE=/home/ron/Bak/dynamic_tunnels/config.txt
rm /home/ron/Bak/dynamic_tunnels/runtime_stats.txt
sed -i '2s/.*/iterations 1000000/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_1Miterations_3SB
mv found_tunnels/* computed_at_night/1TQN_09_1Miterations_3SB
mv runtime_stats.txt computed_at_night/1TQN_09_1Miterations_3SB




