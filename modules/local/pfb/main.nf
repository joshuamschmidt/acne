process MAKEPFB {
    publishDir "$params.outdir/PFB/", pattern: "*.pfb", mode: 'copy'

    label 'pythonTasks'

    input:
    tuple val(id), path(gs_file)

    output:
    tuple val(id), path("${gs_file.baseName}.pfb") , emit: pfb

    script:
    """
    pennCNVtools.py pfb --input $gs_file \
    --output "${gs_file.baseName}.pfb"
    """
}
