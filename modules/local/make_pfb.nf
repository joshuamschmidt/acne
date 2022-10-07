process MAKEPFB {
    publishDir "$params.outdir/PFB/", pattern: "*.pfb", mode: 'copy'

    container 'joshmschmidt/penncnvtools:0.0.1'

    input:
    tuple val(meta), path(gsfile)

    output:
    tuple val(meta), path("${gsfile.baseName}.pfb")

    script:
    """
    pennCNVtools.py pfb --input $gs_file \
    --output "${gsfile.baseName}.pfb"
    """
}
