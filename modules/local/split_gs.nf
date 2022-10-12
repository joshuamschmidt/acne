process SPLITGS {
    container 'joshmschmidt/penncnvtools:0.0.1'

    input:
    tuple val(meta), path(gsfile)

    output:
    tuple val(meta), path("*.txt"), emit: output

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    
    """
    export POLARS_MAX_THREADS=${task.cpus}
    pennCNVtools.py split --input $gsfile
    """
}
