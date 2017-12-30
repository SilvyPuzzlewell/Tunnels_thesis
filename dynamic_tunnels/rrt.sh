#!/bin/bash
if [ -d ./found_tunnels/tunnels0 ]; then
	rm -r ./found_tunnels/tunnels0
fi
mkdir ./found_tunnels/tunnels0

if [ -d ./clusters/clusters0 ]; then
	rm -r ./clusters/clusters0
fi
mkdir ./found_tunnels/tunnels0
mkdir ./clusters/clusters0

cd ./found_tunnels/tunnels0
date > date_pre.txt
cd ~/bak_repository/dynamic_tunnels

./rrt
cd ./found_tunnels/tunnels0
date > date_post.txt
cd ~/Bak/dynamic_tunnels
