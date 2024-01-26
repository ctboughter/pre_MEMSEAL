#!/bin/bash
#$ -N D3001.cite_ranger
#$ -pe threaded 16
cd /hpcdata/lisbcb/MEMSEAL/raw/pool5_data/yfv
module load cellranger/7.1.0
cellranger multi --id=multi_d3001 \
--csv=Dtemp_lib01.csv \
--localcores=16 --localmem=128 > vdj_cited3001.out
