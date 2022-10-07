process SPLITGS {
    container 'joshmschmidt/penncnvtools:0.0.1'

    input:
    tuple val(meta), path(gsfile)

    output:
    tuple val(meta), path("*"), emit: output

    script:
    """
    pennCNVtools.py split --input $gs_file
    """
}
