#!/bin/bash

# So you can either list numbers
# for i in 01 09 14 15
# or count the numbers off
# for i in {01..15} which IS smart enough to go 01 02 .. 14 15
dis=yfv
day=5
# The prefix will change depending on the dataset
# and the assay being moved... But hopefully these 
# general pattterns should be consistent...
prefix=GER
for i in {01..16};
do
    # MAKE SURE THERE ARE NO SPACES
    # IN YOUR VARIABLE DECLARATION
    datdir=gex_"$dis"_run$day$i
    mkdir $datdir

    # Check that _L000 is in your new files...
    # There has been little consistency
    # in the data that we receive.
    mv *"$prefix"GXD$day$i*_L00*/* $datdir/.
    
    # Do NOT remove the files until you know
    # for sure that the script is functioning.
    #rm -r *"$prefix"GXD$day$i*_L00*

done
