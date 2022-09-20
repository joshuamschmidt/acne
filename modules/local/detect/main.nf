process PENNCNV_GC {

    container 'genomicslab/penncnv:1.0.5'

    input:
    tuple val(subBatch), path(pfb)
    path(gc_file)

    output:
    tuple val(subBatch), path("${pfb.baseName}.gc")

    script:

    """
    /root/PennCNV/cal_gc_snp.pl \
    --output "${pfb.baseName}.gc" \
    $gc_file \
    $pfb;
    """
}



/*
process PENNCNV_DETECT {

    container 'genomicslab/penncnv:1.0.5'

    publishDir "$params.outdir/logs", pattern: "*.rawcnv.log", mode: 'copy'
    publishDir "$params.outdir/raw_calls", pattern: "*.rawcnv", mode: 'copy'

    input:
    tuple val(id), path(input), path(pfb), path(hmm), path(gc_model)

    output:
    path "${input.baseName}.rawcnv", emit: cnv
    path "${input.baseName}.rawcnv.log", emit: cnv_log
    
    script:

    """
    detect_cnv.pl \
    -test \
    --loh \
    --confidence \
    --hmmfile $hmm \
    --pfbfile $pfb \
    --logfile ${input.baseName}.rawcnv.log
    --gcmodelfile $gc_model \
    --output  ${input.baseName}.rawcnv \
    $input;
    """

}
*/
