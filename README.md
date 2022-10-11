# acne
nf pipeline for CNV analysis of genotype array data

### 





### Example analysis
Using....

```
nextflow run ./main.nf -resume -profile garvan --input data/anzrag_sample_sheet.txt --partition true --partition_n 300 --hmm /share/ScratchGeneral/jossch/array_cnv/acne/assets/hh550.hmm --gc_model /share/ScratchGeneral/jossch/array_cnv/acne/assets/hg19.gc5Base.txt
```