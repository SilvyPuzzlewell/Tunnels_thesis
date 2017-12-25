#$1 = stats file name
#$2 = save directory

basedir="./relevant/static/$2"

mkdir -p "$basedir/found_tunnels"

python to_latex.py $1 >> latex_stats_code.txt
mv latex_stats_code.txt "$basedir"
cp stats.txt "$basedir"
cp $1 "$basedir"   #stats
cp config.txt "$basedir"
cp -r found_tunnels/* "$basedir/found_tunnels"


