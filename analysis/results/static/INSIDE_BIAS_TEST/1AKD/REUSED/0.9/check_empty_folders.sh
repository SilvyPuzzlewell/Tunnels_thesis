#!/bin/bash
cd ./found_tunnels
for i in {0..100..1}
  do
    if [ -z "$(ls -A tunnels$i)" ]; then
       echo "Empty $i"
    else
       echo "Not Empty"
    fi
 done
cd ..
