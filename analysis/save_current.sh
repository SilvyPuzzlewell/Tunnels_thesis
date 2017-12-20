#$1 = stats file name
#$2 = save directory

basedir="./relevant/static"

mkdir -p "$basedir/$2/found_tunnels"

cp $1 "$basedir/$2"   #stats
cp config.txt "$basedir/$2"
cp -r found_tunnels/* "$basedir/$2/found_tunnels"
