# Processing 'Omics Data for the MEMSEAL DARPA Project

# Summary
This repository should be used as a general reference for the analysis pipeline used for the MEMSEAL DARPA project. All scripts have already been uploaded, and ideally example files can be included to show the "before and after" of running these scripts. Unfortunately, in many cases the files are much too large to host on GitHub.

# Interpreting the Scripts
Scripts are numbered roughly in the order they are run in the pipeline. If there are alternative or parallel processes at a given step, a letter is included (as in step3A, step3B). The entire analysis in this project includes various 'Omics inputs, so eventually this repository will be subdivided for each pipeline. Ordered scripts will be within each directory. For instance, directories will be named "RNAseq", "ATACseq", etc. and each of these directories will have numbered steps.

# Instructions for Use - Bash Scripts for Moving Data and Running Cellranger (Steps 1 to 3)
(Contact Chris if you have questions)

1. The 1move_files.sh script should be run from whichever directory the files are located. Ideally, a given directory will hold data across multiple days, multiple samples, and multiple assays. This script will allow you to move files from whatever poorly formatted directories we get from a sequencing core into a format that is readable by the scripts later on in this pipeline

2. The 2strip_files.sh script is a (potentially optional) way of removing leading stings from the filenames within the directories you created and filled in step 1. Currently it isn't automated to run through all the files. You first run the script with:

```
./2strip_files.sh
```

Then hit enter, and *then* enter the directory name where you'd like to strip the leading strings from the files in said directory.

3. If you don't like how I used these scripts to move/rename files, there's no reason you need to follow this convention. You can follow whatever convention you prefer, as long as you preserve the raw data and end up with H5 files as output (after this step). I've included two bash scripts (you'll run one of them) to automate the process of submitting jobs to a (PBS/SGE) cluster to run cellranger. Basically automating a submission of 16 separate jobs to the cluster. If you won't be running these on a cluster, you can look in the 3outputs directory to see what these scripts output, and just run that locally while changing the relevant parts to run on each individual sample directory.

There are 3 files you need to properly run this script:
- The CITEseq panel used for these experiments (3inputs/citeSeq_codes.csv)
- A Reference Genome (should be on whatever machine you're using)
- A VDJ Reference (also should have already on a cluster)

# Instructions for Use - R Scripts for Normalization and Demultiplexing (Step 4)
(Contact Budha if you have questions)
- The input h5 should be taken from the "outs" directory that you generated in step 3. The names of the files and labels_hashtag excel sheet first column should be identical.

- In the R script the h5 files are read, demultiplexed, SCtransformed and saved on the disk in a loop.

- The column E (fifth column) of the labels_hashtags sheet contains the hashtag information for each index. Notice that the contents of these rows are not identical. They are not identical for two reasons. 1. Experimental design incorporating different hashtags were different on different time points. 2. Some hashtags did not work on some indexes, in an apparently random fashion. These hashes were manually removed.

- How did we know which hashtag did not work for a given label? This is where some manual trial and errors are neccessary. When a given hashtag has failed in a given index, while running the script on that index an error like the following will be seen:

```
Error in HTODemux(D, assay = "HTO", positive.quantile = 0.999, verbose = TRUE) : Cells with zero counts exist as a cluster.
```

- When you encounter such error(s), go back to the code and calculate the total number of reads in each hashtag in a given index with a code something like the following:

```
sums <- rowSums(as.matrix(D@assays$HTO@counts))
```

Remove the hashtag with the minimum reads from the excel sheet (column 5) and run the script again. 
Repeat this for other hashes, with and/or without replacing. The hash with the least reads
is generally, but not always the problematic hashtag. You may have to look at overall distribution
of reads of individial hastags in each index. In extreme scenarios the hashtag quality may be 
really bad and the data may not be usable. In such extereme scenarios we shall meet, go through the 
data and then take a call.

If all goes well you will see something like the following message for a given index:

Cutoff for Ms.Hashtag-1 : 6 reads
Cutoff for Ms.Hashtag-2 : 6 reads
Cutoff for Ms.Hashtag-3 : 5 reads
Cutoff for Ms.Hashtag-4 : 3 reads
Cutoff for Ms.Hashtag-5 : 12 reads
Cutoff for Ms.Hashtag-6 : 3 reads
Cutoff for Ms.Hashtag-7 : 110 reads
Cutoff for Ms.Hashtag-8 : 4 reads
Cutoff for Ms.Hashtag-9 : 6 reads
Cutoff for Ms.Hashtag-10 : 11 reads
Cutoff for Ms.Hashtag-11 : 4 reads
Cutoff for Ms.Hashtag-12 : 26 reads

# Further Reading
The following links may be helpful for understanding the various steps of the analysis:

- Cellranger Multi Instructions: https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-5p-multi
- Seurat Vignettes: https://satijalab.org/seurat/articles/get_started.html
- HTO Demultiplex Docs: https://satijalab.org/seurat/articles/hashing_vignette.html 

# Repository Maintainers
At present, this repository is jointly maintained by Chris Boughter & Budha Chatterjee. Don't hesitate to reach out if there are questions. This project is funded by DARPA, and supported by researchers at the NIH and at the UMD Medical School.