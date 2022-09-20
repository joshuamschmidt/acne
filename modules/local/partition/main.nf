process PARTITIONGS {

    container 'joshmschmidt/penncnvtools:0.0.1'

    input:
    path(gs_file)
    val(partition_size)

    output:
    path("*.partition") , emit: gs

    script:
    """
    pennCNVtools.py partition --input $gs_file --n $partition_size
    """
}
