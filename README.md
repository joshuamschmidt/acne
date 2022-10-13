# acne
nf pipeline for CNV analysis of genotype array data

### 





### Example analysis
Using....

```
export NXF_OPTS="-Xms4G -Xmx8G -Dnxf.pool.maxThreads=2000"

nextflow run ./main.nf -resume -profile garvan \
--input data/anzrag_sample_sheet.txt --partition true \
--partition_n 300 \
--hmm /share/ScratchGeneral/jossch/array_cnv/acne/assets/hh550.hmm \
--gc_model /share/ScratchGeneral/jossch/array_cnv/acne/assets/hg19.gc5Base.txt
```



### TODO
  

Find that there are many NaN values in PFB files - presumably this is because high number of missing genotypes? the GS python scripts
should probably take of care of this, and filter to 95% or 98% genotyping rate. And create a report with the filtered SNPs?

The PENN CNV detect process is called per ind sample. Each of these jobs takes ~15-30secs, though some can take upto 2 mins. Given the cost/overheads of spinning up a compute job on HPC environments, this is probably not optimal. Rather, sub batching into groups of 5-10 samples could probably make sense?

Ind sample outputs do not have batch, job ID etc identifiers appended as prefix/suffix nor are stored in batch/id specific subfolders.  
This should change! have now added prefix to ind split files.
Will need to update per sample meta data to make it ind sample specific.