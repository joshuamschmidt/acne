process PFBTOBED {
    tag "$meta.id"
    label 'process_single'
    container 'joshmschmidt/penncnvtools:0.0.1'


    input:
    tuple val(meta), path(pfb)
    val(gc_window)
    path(chr_sizes)

    output:
    tuple val(meta), path("${pfb.baseName}.bed"), emit: output

    script:
    """
    pfb_to_bed.py \
    --pfb $pfb \
    --window $gc_window \
    --chr_sizes $chr_sizes \
    --output  "${pfb.baseName}.bed"
    """
}
