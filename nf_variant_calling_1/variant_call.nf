params {
    fastqFile = 'path/to/raw.fastq.gz'
    referenceGenome = 'path/to/reference.fasta'
    bwaIndex = 'path/to/reference.fasta'
    gatkPath = 'path/to/gatk.jar'
}

process bwaAlign {
    input:
    file(fastqFile) from fastqChannel

    output:
    file("${fastqFile.baseName}.bam") into bamChannel

    """
    bwa mem -t 4 ${bwaIndex} ${fastqFile} | samtools sort -@ 4 -o ${fastqFile.baseName}.bam
    """
}

process markDuplicates {
    input:
    file bamFile from bamChannel

    output:
    file("${bamFile.baseName}.dedup.bam") into dedupBamChannel

    """
    gatk MarkDuplicatesSpark -I ${bamFile} -O ${bamFile.baseName}.dedup.bam
    """
}

process baseRecalibration {
    input:
    file dedupBamFile from dedupBamChannel

    output:
    file("${dedupBamFile.baseName}.recal.bam") into recalBamChannel

    """
    gatk BaseRecalibrator -I ${dedupBamFile} -R ${referenceGenome} -O ${dedupBamFile.baseName}.recal.table
    gatk ApplyBQSR -I ${dedupBamFile} -R ${referenceGenome} --bqsr-recal-file ${dedupBamFile.baseName}.recal.table -O ${dedupBamFile.baseName}.recal.bam
    """
}

process variantCalling {
    input:
    file recalBamFile from recalBamChannel

    output:
    file("${recalBamFile.baseName}.vcf") into vcfChannel

    """
    gatk HaplotypeCaller -I ${recalBamFile} -R ${referenceGenome} -O ${recalBamFile.baseName}.vcf
    """
}

workflow {
    fastqChannel = Channel.fromPath(params.fastqFile)
    bamChannel = Channel.create()
    dedupBamChannel = Channel.create()
    recalBamChannel = Channel.create()
    vcfChannel = Channel.create()

    bwaAlign()
    markDuplicates()
    baseRecalibration()
    variantCalling()
}
