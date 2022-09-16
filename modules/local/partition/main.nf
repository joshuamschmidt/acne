process PARTITIONGS {

    label 'pythonTasks'

    input:
    tuple val(id), val(partition_size), path(gs_file)

    output:
    tuple val(id), path("*-partition.txt") , emit: gs

    script:
    """
    pennCNVtools.py partition --input $gs_file --n $partition_size
    """
}
