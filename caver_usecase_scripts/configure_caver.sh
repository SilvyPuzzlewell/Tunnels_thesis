#!/bin/bash

rm	~/Bak/caver/examples/QUICK_START/md_snapshots/*
cp	./jd/$1	~/Bak/caver/examples/QUICK_START/md_snapshots/$1
cd	~/Bak/caver/examples/QUICK_START/inputs
./caver.bat
cd ~/bak_repository/caver_usecase_scripts/dynamic_tunnels
