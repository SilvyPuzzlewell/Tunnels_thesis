

#!/bin/bash
#rm -r ./caver_tunnels/*
#python make_atom_pdb.py
python compute_stats_caver_method.py
python process_stats.py $1
