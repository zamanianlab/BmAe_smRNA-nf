#!/usr/bin/env nextflow

// Nextflow.configuration
aux=config.aux_location
data=config.home_location // or brc_location
output=config.output_location
aedesgenome=config.aedesgenome_location

large_core=config.large_core
small_core=config.small_core

// Parameters

params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"


////////////////////////////////////////////////
// ** - Pull in fq files (paired)
////////////////////////////////////////////////

fqs = Channel.fromPath(data + "${params.dir}/*.f[a-z]*q.gz")
                        .map { n -> [ n.getName(), n ] }


////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trim_reads {

   cpus large_core
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.html'
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.json'

   input:
       tuple val(id), file(reads) from fqs

   output:
       tuple id_out, file("${id_out}.fq.gz") into trimmed_fqs
       tuple file("*.html"), file("*.json")  into trim_log

  script:
      id_out = id.replace('.fastq.gz', '')

   """
       fastp -i $reads -o ${id_out}.fq.gz -y -l 15 -h ${id_out}.html -j ${id_out}.json
   """
}
trimmed_fqs.set { trimmed_reads_bwa }


////////////////////////////////////////////////
// ** - LOAD in Aedes genome, indices, and miRNA files
////////////////////////////////////////////////
// aedesgenome points to /mnt/genomes/Other/Aedes_aegypti/

geneset_gtf = file("${aedesgenome}/annotation/geneset_h.gtf.gz")
genome_fa = file("${aedesgenome}/genome.fa")
bwa_indices = Channel.fromPath("${aedesgenome}/BWAIndex/*").collect()


////////////////////////////////////////////////
// ** - bwa mapping
////////////////////////////////////////////////

process align {
    publishDir "${output}/bwa_stats/", mode: 'copy'

    cpus large_core
    tag { id }

    input:
        tuple val(id), file(reads) from trimmed_reads_bwa
        file bwa_indices from bwa_indices.first()

    output:
        file("${id}_align.txt") into bwa_stats

    script:
        fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
        index_base = bwa_indices[0].toString() - ~/.fa*/

        """
        bwa aln -o 0 -n 0 -t ${large_core} ${index_base}.fa ${reads} > ${id}.sai
        bwa samse parasite.fa ${id}.sai ${reads} > ${id}.sam
        samtools view -bS ${id}.sam > ${id}.unsorted.bam
        rm *.sam
        samtools flagstat ${id}.unsorted.bam
        samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
        rm *.unsorted.bam
        samtools index -b ${id}.bam
        samtools flagstat ${id}.bam > ${id}_align.txt
        """
}
