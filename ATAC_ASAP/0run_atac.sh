# Comment crap out to just run as script
!/bin/bash
#$ -N asap1_1
#$ -l mem_free=128G
#$ -l h_vmem=128G

# This is just your basic cellranger script to go from FASTQs into a count matrix
# Only run this if you have already run the "move file" and strip file" scripts

cd /hpcdata/lisbcb/MEMSEAL/ATAC/YFV_atac
module load cellranger-atac/2.0.0

day=1
s=01

cellranger-atac count --id=out_y$day$s \
    --reference ../../ref/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
    --fastqs=atac_yfv_$day$s --sample=ATACYD$day$s > atac$day$s.out
