#!/usr/bin/env nextflow

Channel
    .fromFilePairs( params.reads, size: params.singleend ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}." }
    .set { raw_reads_bbduk } // .into for two or more target channels


/*
 * Step 1. Remove adapters, contaminants and filter low quality
 * bases/reads using bbduk.
 */
process bbduk {
    container 'metashot/bbtools:38.79'
    cpus 2
    memory '4 GB'
    tag "${id}"

    //publishDir "${params.outdir}/bbduk" , mode: 'copy'
    
    input:
    set val(id), path(reads) from raw_reads_bbduk

    output:
    set val(id), path("clean_{1,2}.fastq.gz") into clean_reads_spades

    script:
    """
    bbduk.sh \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out1=noadapt_1.fastq.gz \
        out2=noadapt_2.fastq.gz \
        ref=adapters \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        tpe \
        tbo \
        interleaved=f \
        threads=${task.cpus}
 
    bbduk.sh \
        in1=noadapt_1.fastq.gz \
        in2=noadapt_2.fastq.gz \
        out1=nocont_1.fastq.gz \
        out2=nocont_2.fastq.gz \
        ref=artifacts,phix \
        k=31 \
        hdist=1 \
        interleaved=f \
        threads=${task.cpus}
    
    bbduk.sh \
        in1=nocont_1.fastq.gz \
        in2=nocont_2.fastq.gz \
        out1=clean_1.fastq.gz \
        out2=clean_2.fastq.gz \
        minavgquality=3 \
        maxns=4 \
        qtrim=r \
        trimq=6 \
        mlf=0.5 \
        minlength=${params.minlength} \
        interleaved=f \
        threads=${task.cpus}
    """
}

/*
 * Step 2. Assembly with Spades.
 */
process spades {
    container 'metashot/spades:3.14.0'
    cpus 4
    memory '8 GB'
    tag "${id}"

    publishDir "${params.outdir}/spades" , mode: 'copy',
        saveAs: {filename -> if (filename.indexOf("scaffolds.fasta") > 0) "${id}_scaffolds.fasta"}  

    input:
    set val(id), path(reads) from clean_reads_spades

    output:
    path "spades/*_scaffolds.fasta"

    script:
    task_memory_GB = task.memory.toGiga()
    """
    spades.py \
        --meta \
        --only-assembler \
        -1 ${reads[0]} \
        -2 ${reads[1]}  \
        -k 21,33,55,77,99 \
        --threads ${task.cpus} \
        --memory ${task_memory_GB} \
        -o spades
    """
}