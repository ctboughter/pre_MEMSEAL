#!/bin/bash

# Alright so in this we're going to try to automate submitting these
# Cellranger jobs, instead of the low-throughput I've been doing

for i in {01..16};
do
    echo '#!/bin/bash' > Dtemp_RunRanger$i.sh
    echo "#$ -N D30$i.cite_ranger" >> Dtemp_RunRanger$i.sh
    echo "#$ -pe threaded 16" >> Dtemp_RunRanger$i.sh
    echo "cd /hpcdata/lisbcb/MEMSEAL/raw/pool5_data/yfv" >> Dtemp_RunRanger$i.sh
    echo "module load cellranger/7.1.0" >> Dtemp_RunRanger$i.sh
    echo "cellranger multi --id=multi_d30$i \\" >> Dtemp_RunRanger$i.sh
    echo "--csv=Dtemp_lib$i.csv \\" >> Dtemp_RunRanger$i.sh
    echo "--localcores=16 --localmem=128 > vdj_cited30$i.out" >> Dtemp_RunRanger$i.sh

    head -9 vdj1.csv > Dtemp_lib$i.csv
    echo "GEYGXD30$i,/hpcdata/lisbcb/MEMSEAL/raw/pool5_data/yfv/gex_yfv_runD30$i,1|2|3|4,gene expression," >> Dtemp_lib$i.csv
    echo "CSYGXD30$i,/hpcdata/lisbcb/MEMSEAL/raw/pool5_data/yfv/gex_yfv_runD30$i,1|2|3|4,antibody capture," >> Dtemp_lib$i.csv
    echo "VDYGXD30$i,/hpcdata/lisbcb/MEMSEAL/raw/pool5_data/yfv/vdj_yfv_run30$i,1|2|3|4,vdj," >> Dtemp_lib$i.csv

    qsub Dtemp_RunRanger$i.sh
done
