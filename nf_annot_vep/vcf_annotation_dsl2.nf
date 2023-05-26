#!/usr/bin/env nextflow
nextflow.enable.dsl=2


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
        nextflow run vcf_annotation_dsl2.nf -profile standard --vcf=sample.vcf.gz
        nextflow run vcf_annotation_dsl2.nf -profile docker --vcf=sample.vcf.gz  
        nextflow run vcf_annotation_dsl2.nf -profile batch --vcf=sample.vcf.gz
                 NOTE: use -resume option if nothing changed

    """.stripIndent()
}


process SPLIT_TO_HEADER_AND_BODY{
    tag " Split ${vcf_full} into two separate files: header and body ... "  
    publishDir "${params.outdir}", mode: 'symlink'
    // fing a container with bcftools, now it is local / biocontainers/bcftools:v1.9-1-deb_cv1 
    //container = "broadinstitute/gatk"   
    //containerOptions "-v ${params.inputdir}:/usr/working"   

    input:
        file vcf_full  

    output:
        path "header.vcf", emit: header
        path "no_header.vcf", emit: body

    script:
    """
    bcftools view --header-only ${vcf_full} --output header.vcf
    bcftools view --no-header ${vcf_full} --output no_header.vcf

    """

}


process RUN_VEP {
    tag " VEP annotation of chunk ${vcf_chunk} ... "  
    publishDir "${params.outdir}", mode: 'symlink'
    container = "ensemblorg/ensembl-vep"   // https://hub.docker.com/r/nfcore/base  can also be used
    containerOptions "-v ${params.inputdir}:/usr/working"   

    input:
        file vcf_chunk

    output:
        file "annotated_${vcf_chunk}"

    script:
    """
    vep -i ${vcf_chunk} -o annotated_${vcf_chunk} --database --format "vcf" --no_stats --vcf
    """
    //other possible options: --assembly GRCh38 / use offline database for speed
}



workflow {

    main:
        // Show help message
        if (params.help) {
            helpMessage()
            exit 0
        }

        // full path to input VCF file
        file = params.inputdir + params.vcf

        println """\
            ========================================================
                        VCF  ANNOTATION IN SMALL CHUNKS  
            ========================================================
            input VCF to annotate     : ${file}
            out folder                : ${params.outdir}
            note: use -resume option if nothing changed

            """
            .stripIndent()

        ch_vcf = Channel.fromPath(file)

        splitted_vcf = SPLIT_TO_HEADER_AND_BODY(ch_vcf)

        ch_chunks = splitted_vcf.body.splitText(by:2000, compress:true, file:true)
        //ch_chunks = Channel.fromPath(file).splitText(by:500, decompress:true, compress:true, file:true)

        ch_annotated_chunks = RUN_VEP(ch_chunks)

    emit:
        // TODO: how keep the order of batches ? + need to bgzip output
        ch_annotated_chunks.
            concat(splitted_vcf.header).
            collectFile(name: "${params.outdir}result.vcf", newLine: true)
 


}
