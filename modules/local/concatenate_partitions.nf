process CONCATENATE_PARTITIONS {
    tag "$meta.id"
    //container 'genomicslab/penncnv:latest'
    //container 'wallen/penncnv:1.0.5'

    publishDir "$params.outdir/final_penn_calls", pattern: "${prefix}.filteredcnv", mode: 'copy'
    publishDir "$params.outdir/final_penn_calls", pattern: "${prefix}.qcpass", mode: 'copy'
    publishDir "$params.outdir/final_penn_calls", pattern: "${prefix}.qcsum", mode: 'copy'

    input:
    tuple val(meta), path('*.filteredcnv'), path('*.qcpass'), path('*.qcsum')

    output:
    tuple val(meta), path ("${prefix}.filteredcnv"), path ("${prefix}.qcpass"), path ("${prefix}.qcsum"), emit: output
    
    script:
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    cat *.filteredcnv > "$prefix".filteredcnv
    cat *.qcpass > "$prefix".qcpass
    cat *.qcsum > "$prefix".qcsum
    """
}