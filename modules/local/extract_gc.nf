process EXCTRACTGC {
    tag "$meta.id"
    label 'process_single'
    container 'joshmschmidt/penncnvtools:0.0.1'


    input:
    tuple val(meta), path(nuc)

    output:
    tuple val(meta), path("${nuc.baseName}.gc"), emit: gc

    script:
    """
    extract_gc.py \
    --nuc $nuc \
    --output  "${nuc.baseName}.gc"
    """
}
