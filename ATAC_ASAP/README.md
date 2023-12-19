# Specific Processing for ATACseq

# Quick Notes
So more comprehensive notes will follow at some point, but for now, mostly posting so Budha can help me with the R outputs.

I recommend downloading from the MEMSEAL directory the entire ATAC_YFV_303 output directory at this path:

```
/hpcdata/lisbcb/MEMSEAL/ATAC/YFV_atac/out_y303/outs
```

Exclude the .bam files in this directory though, because these are 18GB (the largest files by far), and you don't need them for analysis

The script runs fine as-is, assuming you can get signac working, the outputs are just stripped of the annotations for the chromosomal data (i.e. the counts for "Cd3d" lose the "Cd3d" label and revert back to "chr2-87011729-87035519"). This, and I also wasn't able to export the aligned reads for coverage plots (i.e. the "peaks" in the dataframe). This might have something to do with the final conversion step I'm running, I'm not sure. But that step is necessary.

p.s. Hopefully signac just installs alright on your machine. It failed to install on 3 separate machines I'm using... If it does fail and you have issues getting docker/signac working, then let me know.