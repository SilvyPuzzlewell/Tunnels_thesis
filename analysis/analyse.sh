./run_analysis.sh
./save_current.sh processed_stats.txt $1
./run_analysis_norm.sh

basedir="./results/static/$1"
mv processed_stats.txt "$basedir/processed_stats_norm.txt"
mv iterations_tunnel_counts.txt "$basedir/iterations_tunnel_counts_norm.txt" 
