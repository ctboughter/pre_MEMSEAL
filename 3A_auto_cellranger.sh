#!/bin/bash

# Alright so in this we're going to try to automate submitting these
# Cellranger jobs, instead of the low-throughput I've been doing

d=3

for i in {01..16};
do
    echo '#!/bin/bash' > Dtemp_RunRanger$d$i.sh
    echo "#$ -N D$d.rub$i.cite_ranger" >> Dtemp_RunRanger$d$i.sh
    echo "#$ -pe threaded 16" >> Dtemp_RunRanger$d$i.sh
    echo "cd /hpcdata/lisbcb/MEMSEAL/pool3_data" >> Dtemp_RunRanger$d$i.sh
    echo "module load cellranger/7.1.0" >> Dtemp_RunRanger$d$i.sh
    echo "cellranger count --id=gex_cite_rubD$d$i \\" >> Dtemp_RunRanger$d$i.sh
    echo "--transcriptome=/hpcdata/bio_data/cellranger/refdata-gex-mm10-2020-A \\" >> Dtemp_RunRanger$d$i.sh
    echo "--libraries=Dtemp_lib$d$i.csv \\" >> Dtemp_RunRanger$d$i.sh
    echo "--feature-ref=citeSeq_codes.csv \\" >> Dtemp_RunRanger$d$i.sh
    echo "--localcores=16 --localmem=128 > gex_cited$d$i.out" >> Dtemp_RunRanger$d$i.sh

    echo "fastqs,sample,library_type," > Dtemp_lib$d$i.csv
    echo "/hpcdata/lisbcb/MEMSEAL/pool3_data/gex_rub_runD$d$i,GERGXD$d$i,Gene Expression," >> Dtemp_lib$d$i.csv
    echo "/hpcdata/lisbcb/MEMSEAL/pool3_data/gex_rub_runD$d$i,CSRGXD$d$i,Antibody Capture," >> Dtemp_lib$d$i.csv

    qsub Dtemp_RunRanger$d$i.sh
done
