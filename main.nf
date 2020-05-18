#!/usr/bin/env nextflow

Channel
    .fromFilePairs( params.reads, size: (params.single_end || params.interleaved) ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}." }
    .set { raw_reads } // .into for two or more target channels

/*
 * Step 0. Deinterleave paired reads.
 */
if (params.interleaved) {
    process reformat {
        container 'metashot/bbtools:38.79-1'
        
        cpus 2
        memory { 2.GB * (2**(task.attempt-1)) }
        maxRetries 3
        errorStrategy 'retry'

        tag "${id}"
    
        input:
        tuple val(id), path(reads) from raw_reads

        output:
        tuple val(id), path("read_*.fastq.gz") into raw_reads_stats, raw_reads_adapter_trimming

        script:
        task_memory_GB = task.memory.toGiga()
        // Note: Even with "t=1", reformat will generally use over 2 CPU cores on average 
        // since the I/O is in separate threads.
        """
        reformat.sh \
            -Xmx${task_memory_GB}g \
            in=$reads \
            out1=read_1.fastq.gz \
            out2=read_2.fastq.gz \
            t=1
        """
    }
} else {
    raw_reads.into { raw_reads_stats; raw_reads_adapter_trimming }
}

/*
 * Step 1. Raw reads histograms.
 */

process raw_reads_stats {
    container 'metashot/bbtools:38.79-1'
    
    cpus 2
    memory { 2.GB * (2**(task.attempt-1)) }
    maxRetries 3
    errorStrategy 'retry'

    tag "${id}"

    publishDir "${params.outdir}/samples/${id}/raw_reads_stats" , mode: 'copy'

    input:
    tuple val(id), path(reads) from raw_reads_stats

    output:
    path "*hist.txt"

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    """
    bbduk.sh \
        -Xmx${task_memory_GB}g \
        $input \
        bhist=bhist.txt \
        qhist=qhist.txt \
        gchist=gchist.txt \
        aqhist=aqhist.txt \
        lhist=lhist.txt \
        gcbins=auto
    """
}
 
/*
 * Step 2.a Remove adapters.
 */
if (!params.skip_adapter_trimming) {
    process adapter_trimming {
        container 'metashot/bbtools:38.79-1'
        
        cpus 4
        memory { 2.GB * (2**(task.attempt-1)) }
        maxRetries 3
        errorStrategy 'retry'

        tag "${id}"
    
        publishDir "${params.outdir}/samples/${id}/bbduk" , mode: 'copy',
            pattern: "stats_adapter_trimming.txt"

        input:
        tuple val(id), path(reads) from raw_reads_adapter_trimming

        output:
        tuple val(id), path("clean_adapter*.fastq.gz") into clean_adapter_reads
        path "stats_adapter_trimming.txt"

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean_adapter.fastq.gz" : "out1=clean_adapter_1.fastq.gz out2=clean_adapter_2.fastq.gz"
        """
        bbduk.sh \
            -Xmx${task_memory_GB}g \
            $input \
            $output \
            ref=adapters \
            ktrim=r \
            k=23 \
            mink=11 \
            hdist=1 \
            tpe \
            tbo \
            interleaved=f \
            stats=stats_adapter_trimming.txt \
            threads=${task.cpus}
        """
    }
} else {
    raw_reads_adapter_trimming.set { clean_adapter_reads }
}

/*
 * Step 2.b Remove contaminants.
 */
if (!params.skip_contaminant_filtering) {
    process contaminant_filtering {
        container 'metashot/bbtools:38.79-1'

        cpus 4
        memory { 2.GB * (2**(task.attempt-1)) }
        maxRetries 3
        errorStrategy 'retry'

        tag "${id}"
    
        publishDir "${params.outdir}/samples/${id}/bbduk" , mode: 'copy', pattern: "stats_contaminant_filtering.txt"

        input:
        tuple val(id), path(reads) from clean_adapter_reads

        output:
        tuple val(id), path("clean_contaminant*.fastq.gz") into clean_contaminant_reads
        path "stats_contaminant_filtering.txt"

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean_contaminant.fastq.gz" : "out1=clean_contaminant_1.fastq.gz out2=clean_contaminant_2.fastq.gz"
        """
        bbduk.sh \
            -Xmx${task_memory_GB}g \
            $input \
            $output \
            ref=artifacts,phix \
            k=31 \
            hdist=1 \
            interleaved=f \
            stats=stats_contaminant_filtering.txt \
            threads=${task.cpus}
        """
    }
} else {
    clean_adapter_reads.set { clean_contaminant_reads }
}

/*
 * Step 2.c Quality filtering/trimming and length filtering.
 */
if (!params.skip_quality_trimming) {
    process quality_trimming {
        container 'metashot/bbtools:38.79-1'
        
        cpus 4
        memory { 2.GB * (2**(task.attempt-1)) }
        maxRetries 3
        errorStrategy 'retry'

        tag "${id}"

        input:
        tuple val(id), path(reads) from clean_contaminant_reads

        output:
        tuple val(id), path("clean*.fastq.gz") into clean_reads_stats, clean_reads_spades, clean_reads_megahit
  
        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean.fastq.gz" : "out1=clean_1.fastq.gz out2=clean_2.fastq.gz"
        """
        bbduk.sh \
            $input \
            $output \
            maq=10 \
            maxns=4 \
            qtrim=r \
            trimq=6 \
            mlf=0.5 \
            minlen=${params.min_read_len} \
            interleaved=f \
            threads=${task.cpus}
            """
    }
} else {
    clean_contaminant_reads.into { clean_reads_stats; clean_reads_spades; clean_reads_megahit }
}


/*
 * Step 3. Clean reads histograms.
 */

process clean_reads_stats {
    container 'metashot/bbtools:38.79-1'
    
    cpus 2
    memory { 2.GB * (2**(task.attempt-1)) }
    maxRetries 3
    errorStrategy 'retry'

    tag "${id}"

    publishDir "${params.outdir}/samples/${id}/clean_reads_stats" , mode: 'copy'

    when:
    ! (params.skip_adapter_trimming && 
       params.skip_contaminant_filtering &&
       params.skip_quality_trimming)

    input:
    tuple val(id), path(reads) from clean_reads_stats

    output:
    path "*hist.txt"

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    """
    bbduk.sh \
        -Xmx${task_memory_GB}g \
        $input \
        bhist=bhist.txt \
        qhist=qhist.txt \
        gchist=gchist.txt \
        aqhist=aqhist.txt \
        lhist=lhist.txt \
        gcbins=auto
    """
}

/*
 * Step 4.a Assembly with Spades.
 */
if (!params.single_end && !params.megahit_only) {
    process spades {
        container 'metashot/spades:3.14.0-1'
        
        cpus 4
        time = 12.h
        memory { 4.GB * (2**(task.attempt-1)) }
        maxRetries 1
        errorStrategy 'retry'
    
        tag "${id}"
    
        publishDir "${params.outdir}/samples/${id}" , mode: 'copy'
        publishDir "${params.outdir}/samples/${id}" , mode: 'copy' ,
            saveAs: {filename -> if (filename == "spades/scaffolds.fasta") "assembly/scaffolds.fasta"}
    
        input:
        tuple val(id), path(reads) from clean_reads_spades
    
        output:
        tuple val(id), path("spades/scaffolds.fasta") into scaffolds_spades
        path "spades/contigs.fasta"
        path "spades/spades.log"
    
        script:
        task_memory_GB = task.memory.toGiga()
        """
        spades.py \
            --meta \
            --only-assembler \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            --threads ${task.cpus} \
            --memory ${task_memory_GB} \
            -o spades
        """
    }
} else {
    scaffolds_spades = Channel.empty()
}

/*
 * Step 4.b Assembly with Megahit.
 */

if (params.single_end || params.megahit_only) {
    process megahit {
        container 'metashot/megahit:1.2.9-1'

        cpus 4
        time = 12.h
        memory { 1.GB * (2**(task.attempt-1)) }
        maxRetries 3
        errorStrategy 'retry'

        tag "${id}"

        publishDir "${params.outdir}/samples/${id}" , mode: 'copy'
        publishDir "${params.outdir}/samples/${id}" , mode: 'copy' ,
            saveAs: {filename -> if (filename == "megahit/final.contigs.fa") "assembly/scaffolds.fasta"}

        input:
        tuple val(id), path(reads) from clean_reads_megahit

        output:
        tuple val(id), path("megahit/final.contigs.fa") into scaffolds_megahit

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "-r \"$reads\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
        """
        megahit \
            $input \
            -t ${task.cpus} \
            --k-min 27 \
            --k-max 99 \
            --k-step 14 \
            --kmin-1pass \
            --memory $task_memory_GB \
            -o megahit
        """
    }
} else {
    scaffolds_megahit = Channel.empty()
}

scaffolds_spades
    .mix(scaffolds_megahit)
    .set {scaffolds_stats}

/*
 * Step 5. Scaffold statistics.
 */
process assembly_stats {
    container 'metashot/bbtools:38.79-1'
    
    // Stats is single threaded. Stats uses 120MB of RAM regardless
    // of the assembly size.
    cpus 1
    memory 1.GB
 
    tag "${id}"

    publishDir "${params.outdir}/samples/${id}/assembly" , mode: 'copy'

    input:
    tuple val(id), path(scaffolds) from scaffolds_stats

    output:
    path "stats.txt"

    script:
    """
    stats.sh \
        in=$scaffolds \
        out=stats.txt
    """
}
