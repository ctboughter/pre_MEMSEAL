# Add header if need be to run on cluster:
#$ -N proc_asap3_3
#$ -l mem_free=128G
#$ -l h_vmem=128G

# Alright either way we need to run this script:
# We are calling an R script from a singularity container,
# if that all makes sense...

cd /hpcdata/lisbcb/MEMSEAL/ATAC/YFV_atac

# Can use sed to automate running the script
# Should probably make a dummy file with XXX instead of the number
# and then create a new file based on this dummy, so we can run all of
# these at the exact same time.

#sed -i "s/303/304/g" process_ATAC.R

singularity exec --bind /hpcdata/lisbcb/MEMSEAL/ATAC/YFV_atac signac_latest.sif Rscript process_ATAC.R > proc_asap303.out
