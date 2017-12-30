

rm -r ~/bak_repository/analysis/found_tunnels/*
rm -r ~/bak_repository/analysis/clusters/*
cp -r ~/bak_repository/dynamic_tunnels/found_tunnels/* ~/bak_repository/analysis/found_tunnels
cp -r ~/bak_repository/dynamic_tunnels/clusters/* ~/bak_repository/analysis/clusters
cp ~/bak_repository/dynamic_tunnels/config.txt ~/bak_repository/analysis/config.txt
