process SPLITGS {
    container 'joshmschmidt/penncnvtools:0.0.1'

    input:
    path(gs_file)

    output:
    tuple val("${gs_file.baseName}"), path("*.txt")

    script:
    """
    pennCNVtools.py split --input $gs_file
    """
}
