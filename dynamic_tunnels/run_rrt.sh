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
date >> date_pre.log
mv date_pre.log "found_tunnels/date_pre.log"
rm time_log.log
for i in {0..100..1}
do
	echo "iteration $i"
	date >> time_log.log 
	./rrt $i
	date >> time_log.log
	mv runtime_stats.txt "found_tunnels/tunnels$i/runtime_stats.txt" 	
done
date >> date_post.log
mv date_post.log "found_tunnels/date_post.log"
mv time_log.log "found_tunnels/time_log.log"
