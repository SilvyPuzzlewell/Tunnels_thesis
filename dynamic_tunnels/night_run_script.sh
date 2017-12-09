



FILE=/home/fif/bak_repository/dynamic_tunnels/config.txt
rm /home/fif/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '11s/.*/inside_sampling_bias 1' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_500kIB1
mv found_tunnels/* computed_at_night/1TQN_09_500kIB1
mv runtime_stats.txt computed_at_night/1TQN_09_500kIB1
cp config.txt computed_at_night/1TQN_09_500kIB1


FILE=/home/fif/bak_repository/dynamic_tunnels/config.txt
rm /home/fif/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '11s/.*/inside_sampling_bias 2/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_500kIB2
mv found_tunnels/* computed_at_night/1TQN_09_500kIB2
mv runtime_stats.txt computed_at_night/1TQN_09_500kIB2
cp config.txt computed_at_night/1TQN_09_500kIB2


FILE=/home/fif/bak_repository/dynamic_tunnels/config.txt
rm /home/fif/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '11s/.*/inside_sampling_bias 4/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_500kIB4
mv found_tunnels/* computed_at_night/1TQN_09_500kIB4
mv runtime_stats.txt computed_at_night/1TQN_09_500kIB4
cp config.txt computed_at_night/1TQN_09_500kIB4

FILE=/home/fif/bak_repository/dynamic_tunnels/config.txt
rm /home/fif/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '11s/.*/inside_sampling_bias 6/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_500kIB6
mv found_tunnels/* computed_at_night/1TQN_09_500kIB6
mv runtime_stats.txt computed_at_night/1TQN_09_500kIB6
cp config.txt computed_at_night/1TQN_09_500kIB6


FILE=/home/fif/bak_repository/dynamic_tunnels/config.txt
rm /home/fif/bak_repository/dynamic_tunnels/runtime_stats.txt
sed -i '11s/.*/inside_sampling_bias 8/' $FILE
./run_rrt.sh

mkdir computed_at_night/1TQN_09_500kIB8
mv found_tunnels/* computed_at_night/1TQN_09_500kIB8
mv runtime_stats.txt computed_at_night/1TQN_09_500kIB8
cp config.txt computed_at_night/1TQN_09_500kIB8




