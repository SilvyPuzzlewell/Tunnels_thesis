./run_analysis.sh
./save_current.sh processed_stats.txt $1
./run_analysis_norm.sh

basedir="./relevant/static/$1"
mv processed_stats.txt "$basedir/processed_stats_norm.txt" 
