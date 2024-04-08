process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path samplesheet

    output:
    path ("${samplesheet.baseName}-valid.csv")       , emit: csv
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    input_sample_sheet_check.pl \
    --samplesheet $samplesheet \
    --outfile ${samplesheet.baseName}-valid.csv
    """
}
        