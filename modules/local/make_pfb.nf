process MAKEPFB {
    tag "$meta.id"
    publishDir "$params.outdir/PFB/", pattern: "*.pfb", mode: 'copy'
    
    container 'joshmschmidt/penncnvtools:0.0.1'


    input:
    tuple val(meta), path(gsfile)

    output:
    tuple val(meta), path("${gsfile.baseName}.pfb"), emit: output

    script:
    """
    export POLARS_MAX_THREADS=${task.cpus}
    pennCNVtools.py pfb --input $gsfile \
    --output "${gsfile.baseName}.pfb"
    """
}
