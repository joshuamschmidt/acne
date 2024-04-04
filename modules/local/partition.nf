process PARTITIONGS {
    tag "$meta.id"
    container 'joshmschmidt/penncnvtools:0.0.1'

    input:
    tuple val(meta), path(gsfile), path(sample_include)
    val(partition_size)

    output:
    tuple val(meta), path("${prefix}*"), emit: output
    //path '*.partition'

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    export POLARS_MAX_THREADS=${task.cpus}
    pennCNVtools.py partition --input $gsfile --n $partition_size
    """
}
