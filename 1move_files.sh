#!/bin/bash

# So you can either list numbers
# for i in 01 09 14 15
# or count the numbers off
# for i in {01..15} which IS smart enough to go 01 02 .. 14 15
for i in {09..12};
do
    # MAKE SURE THERE ARE NO SPACES
    # IN YOUR VARIABLE DECLARATION
    #datdir=vdj_rubella_run5$i
    datdir=gex_rubella_run5$i
    mkdir $datdir

    #mv *VDRGXD5$i*_L00*/* $datdir/.
    #rm -r *VDRGXD5$i*_L00*
    
    mv *GERGXD5$i*_L00*/* $datdir/.
    #mv *CSRGXD30$i*_L00*/* $datdir/.
    rm -r *GERGXD5$i*_L00*
    #rm -r *CSRGXD30$i*_L00*

done
