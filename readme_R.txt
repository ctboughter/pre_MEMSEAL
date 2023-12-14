The h5 files are to be kept in a directory. The names of the files 
and labels_hashtag excel sheet first column should be identical.

In the R script the h5 files are read, demultiplexed, SCtransformed and 
saved on the disk in a loop.

The column E (fifth column) of the labels_hashtags sheet contains the 
hashtag information for each index. Notice that the contents of these rows 
are not identical. They are not identical for two reasons. 1. Experimental
design incorporating different hashtags were different on different time points.
2. Some hashtags did not work on some indexes, in an apparently random fashion. 
These hashes were manually removed.

How did we know which hashtag did not work for a given label? This is where 
some manual trial and errors are neccessary. When a given hashtag has failed in a given index,
while running the script on that index an error like the following will be seen:

Error in HTODemux(D, assay = "HTO", positive.quantile = 0.999, verbose = TRUE) : 
Cells with zero counts exist as a cluster.

When you encounter such error(s), go back to the code and calculate the total number of
reads in each hashtag in a given index with a code something like the following:

sums <- rowSums(as.matrix(D@assays$HTO@counts))

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

