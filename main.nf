#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow             = params.workflow
idmap                = file(params.idmap, tyep:'file')
outdir               = file(params.outdir, tyep:'dir')
outdir.mkdir()

process run_cram2bam {
    tag { "${input.baseName}" }
    publishDir "${outdir}", mode: 'copy', overwrite: true

    input:
    path(input)
    
    output:
    path("${input.baseName}*.{bam,bam.bai}"), emit: bams

    """
    samtools view -@ ${task.cpus} -b -o ${input.baseName}.bam ${input}
    samtools index -@ ${task.cpus} ${input.baseName}.bam
    """
}

process run_cram2bam_ddd {
    tag { "${cram.baseName}" }
    publishDir "${outdir}/bams", mode: 'copy', overwrite: false
    errorStrategy 'ignore'

    input:
    path(cram)
    
    output:
    path("*.{bam,bam.bai}"), emit: bams

    """
    export REF_PATH=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s

    name=`awk '\$1 == \"${cram.baseName}\" { print \$3 }' ${idmap}`
    samtools view -@ ${task.cpus} -b -o tmp.bam ${cram}
    samtools view -@ ${task.cpus} -H tmp.bam | sed -e 's/SM:\\(.*\\)	CN:/SM:'"\$name"'	CN:/' | \
        samtools reheader - tmp.bam > \$name.bam
    samtools index -@ 5 \$name.bam
    rm tmp.bam
    """
}

process run_sort_cram {
    tag { "${input.baseName}" }

    input:
    tuple val(sample), path(input)
    
    output:
    tuple val(sample), path("${sample}_sorted.bam"), emit: sorted_bam

    """
    export REF_PATH=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s

    samtools sort -n -@ ${task.cpus} ${input} -o ${sample}_sorted.bam
    """
}

process run_cram2fastq {
    tag { "${sample}" }
    
    publishDir "${outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(input)
    
    output:
    tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz"), emit: fastq
    tuple val(sample), path("${sample}_singletons.fastq.gz"), emit: singletons
    tuple val(sample), path("${sample}.fastq.gz"), emit: other

    """
    samtools fastq -@ ${task.cpus} \
            -T "*" -O \
            -0 ${sample}.fastq.gz \
            -1 ${sample}_R1.fastq.gz \
            -2 ${sample}_R2.fastq.gz \
            -s ${sample}_singletons.fastq.gz \
            -n ${input}
    """
}

process run_cram2fastq_BQSR {
    tag { "${sample}" }
    publishDir "${outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(input)
    
    output:
    tuple val(sample), path("${sample}_R1_amended.fastq.gz"), path("${sample}_R2_amended.fastq.gz"), emit: fastq
    tuple val(sample), path("${sample}_amended.singletons.fastq.gz"), emit: singletons
    tuple val(sample), path("${sample}_amended.fastq.gz"), emit: other

    """
    export REF_PATH=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s

    samtools collate -@ ${task.cpus} -u -O ${input} |\
        samtools fastq -@ ${task.cpus} \
            -T "*" -O \
            -0 ${sample}_amended.fastq.gz \
            -1 ${sample}_R1_amended.fastq.gz \
            -2 ${sample}_R2_amended.fastq.gz \
            -s ${sample}_amended_singletons.fastq.gz \
            -n
    """
}

workflow CRAM2BAM {
    take:
    input

    main:
    run_cram2fastq(input)
}

workflow CRAM2BAM_DDD {
    take:
    input

    main:
    run_cram2bam_ddd(input)
}

workflow CRAM2FASTQ {
    take:
    input

    main:
    run_sort_cram(input)
    run_cram2fastq(run_sort_cram.out.sorted_bam)
}

workflow CRAM2FASTQ_BQSR {
    take:
    input

    main:
    run_sort_cram(input)
    run_cram2fastq_BQSR(run_sort_cram.out.sorted_bam)
}

workflow {
    switch (workflow) {
        case['cram2bam']:
            Channel.fromPath(params.input + "/*.cram", type:'file').set { input }
            input.view()
            CRAM2BAM(input)
            break
            // =====
        case['cram2bam_ddd']:
            Channel.fromPath(params.input + "/*.cram", type:'file').set { input }
            CRAM2BAM_DDD(input)
            break
            // =====
        case['bam2cram']:
            break
            // =====            
        case['cram2fastq']:
            Channel.fromPath(params.input + "/*.{cram,bam}", type:'file').set { input }
            CRAM2FASTQ(input)
            break
            // =====
        case['cram2fastq_bqsr']:
            Channel.fromPath(params.samplesheet_crams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}",
                               "${row.CRAM}"] }
                .set { input }
            CRAM2FASTQ_BQSR(input)
            break
            // =====
        default:
            exit 1, "NO WORKFLOW GIVEN!"
            break
            // =====
    }
}
