#!/bin/bash
cd ./found_tunnels
for i in {0..100..1}
  do
  	echo $i
	DIRECTORY="tunnels$i"
  	echo $DIRECTORY
    if [ -d "$DIRECTORY" ] 
    then
    	rm -r ./$DIRECTORY/*
    else
    	mkdir ./$DIRECTORY 
    fi
 done
cd ..
date
rm runtime_stats.txt
for i in {0..100..1}
  do
    echo "iteration $i"
  	./rrt $i
    mv runtime_stats.txt "found_tunnels/tunnels$i/runtime_stats.txt" 	
  done
date
