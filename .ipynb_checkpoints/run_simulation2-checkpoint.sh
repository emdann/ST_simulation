#!/bin/bash

seed=$1
n_spots=$2
id=$3

python assemble_composition_2.py ${seed} --tot_spots $n_spots --assemble_id $id
python assemble_st_2.py ${seed} --assemble_id $id
