#!/bin/bash
rm -r ./found_tunnels/tunnels0
mkdir ./found_tunnels/tunnels0

cd ./found_tunnels/tunnels0
date > date_pre.txt
cd ~/Bak/dynamic_tunnels

./rrt
cd ./found_tunnels/tunnels0
date > date_post.txt
cd ~/Bak/dynamic_tunnels
