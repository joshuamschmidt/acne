process PENNCNV_GC {

    //container 'genomicslab/penncnv:1.0.5'
    container 'wallen/penncnv:1.0.5'
    input:
    tuple val(meta), path(pfb)
    path(gcfile)

    output:
    tuple val(meta), path("${pfb.baseName}.gc"), emit: output

    script:

    """
    cal_gc_snp.pl \
    --output "${pfb.baseName}.gc" \
    $gcfile \
    $pfb;
    """
}




