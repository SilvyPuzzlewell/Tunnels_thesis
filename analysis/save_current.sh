#$1 = stats file name
#$2 = save directory

basedir="./results/static/$2"

mkdir -p "$basedir/found_tunnels"

python to_latex.py $1 >> latex_stats_code.txt
mv latex_stats_code.txt "$basedir"
mv runtime.log "$basedir"
mv times.log "$basedir"
mv iterations_tunnel_counts.txt "$basedir"
cp stats.txt "$basedir"
cp $1 "$basedir"   #stats
cp config.txt "$basedir"
cp -r found_tunnels/* "$basedir/found_tunnels"


