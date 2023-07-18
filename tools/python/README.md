`mamba create --name acnePy python polars numpy zstandard`

`mamba activate acnePy`

`pip list --format=freeze > requirements.txt`

`docker build -t joshmschmidt/penncnvtools:0.0.1 .`

`docker push joshmschmidt/penncnvtools:0.0.1`

`singularity build /share/ClusterShare/software/contrib/jossch/singularity_images/penncnvtools:0.0.1 docker://joshmschmidt/penncnvtools:0.0.1`
