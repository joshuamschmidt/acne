# acne
nf pipeline for CNV analysis of genotype array data

### 





### Example analysis
Using....

```
export NXF_OPTS="-Xms24G -Xmx24G -Dnxf.pool.maxThreads=2000"
nextflow run ./main.nf -resume -profile garvan --input data/anzrag_sample_sheet.txt --partition true --partition_n 300 --hmm /share/ScratchGeneral/jossch/array_cnv/acne/assets/hh550.hmm --gc_model /share/ScratchGeneral/jossch/array_cnv/acne/assets/hg19.gc5Base.txt
```



### TODO
python polars is used for all steps involving the raw GS output. Using optimised lazy evaluation,  
duplicate col names are very problematic. If present causes kernal panic and job failure. In the future, will  
add a process to check format of GS file. Ideally, one could also update the col names e.g. "sample_X and sample_X-DUP"  

Find that there are many NaN values in PFB files - presumably this is because high number of missing genotypes?

polars is also a hungry mofo - and will consume more cpus than requested if able to.  
Look into setting relevant env var to stop this.

The PENN CNV detect process is called per ind sample. Each of these jobs takes ~15-30secs, though some can take upto 2 mins. Given the cost/overheads of spinning up a compute job on HPC environments, this is probably not optimal. Rather, sub batching into groups of 5-10 samples could probably make sense?

Ind sample outputs do not have batch, job ID etc identifiers appended as prefix/suffix nor are stored in batch/id specific subfolders.  
This should change!