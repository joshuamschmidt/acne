process PARTITIONGS {

    container 'joshmschmidt/penncnvtools:0.0.1'

    input:
    tuple val(meta), path(gsfile)
    val(partition_size)

    output:
    tuple val(meta), path '*.partition'

    script:
    """
    pennCNVtools.py partition --input $gsfile --n $partition_size
    """
}
