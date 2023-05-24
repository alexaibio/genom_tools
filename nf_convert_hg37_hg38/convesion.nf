params {
    inputVcf = 'path/to/input.vcf'
    chainFile = 'path/to/hg19ToHg38.over.chain.gz'
    liftoverTool = 'path/to/liftOver'
}

process liftover {
    input:
    file(inputVcf) from inputVcfChannel

    output:
    file("${inputVcf.baseName}.hg38.vcf") into outputVcfChannel

    """
    ${liftoverTool} ${inputVcf} ${chainFile} ${inputVcf.baseName}.hg38.unmapped \
        ${inputVcf.baseName}.hg38.vcf
    """
}

workflow {
    inputVcfChannel = Channel.fromPath(params.inputVcf)
    outputVcfChannel = Channel.create()

    liftover()
}
