process MAKEPFB {
    publishDir "$params.outdir/PFB/", pattern: "*.pfb", mode: 'copy'

    container 'joshmschmidt/penncnvtools:0.0.1'

    input:
    path(gs_file)

    output:
    tuple val("${gs_file.baseName}"), path("${gs_file.baseName}.pfb")

    script:
    """
    pennCNVtools.py pfb --input $gs_file \
    --output "${gs_file.baseName}.pfb"
    """
}
