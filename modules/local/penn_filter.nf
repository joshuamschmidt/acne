process PENNCNV_FILTER {

    //container 'genomicslab/penncnv:latest'
    container 'wallen/penncnv:1.0.5'

    publishDir "$params.outdir/merge_calls", pattern: "*.filteredcnv", mode: 'copy'
    publishDir "$params.outdir/penn_qc", pattern: "*.qcpass", mode: 'copy'
    publishDir "$params.outdir/penn_qc", pattern: "*.qcsum", mode: 'copy'

    input:
    tuple val(meta), path(calls), path(logs)

    output:
    tuple val(meta), path ("${prefix}.filteredcnv"), path ("${prefix}.qcpass"), path ("${prefix}.qcsum"), emit: output
    
    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    filter_cnv.pl  \
    $calls \
    -qclogfile $logs \
    --qclrrsd 0.3 \
    --qcbafdrift 0.01 \
    --qcwf 0.05 \
    -qcpassout ${prefix}.qcpass \
    -qcsumout ${prefix}.qcsum \
    -out ${prefix}.filteredcnv;
    """

}