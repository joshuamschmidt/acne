`mamba create --name acnePy python=3.10.6 polars numpy`

`mamba activate acnePy`

`pip list --format=freeze > requirements.txt`

`docker build -t joshmschmidt/penncnvtools:0.0.1 .`

`docker push joshmschmidt/penncnvtools:0.0.1`

`singularity build /share/ClusterShare/software/contrib/jossch/singularity_images/illumina_gs_to_cnv:0.0.1  docker://joshmschmidt/illumina_gs_to_cnv:0.0.1`
