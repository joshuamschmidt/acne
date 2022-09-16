process SPLITGS {

    label 'pythonTasks'

    input:
    tuple val(id), path(gs_file)

    output:
    tuple val(id), path("*.txt") , emit: pfb

    script:
    """
    pennCNVtools.py split --input $gs_file
    """
}
