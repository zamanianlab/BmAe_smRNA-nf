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

fqs = Channel.fromPath(data + "${params.dir}/*.fastq.gz")
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
       set val(id), file(reads) from fqs

   output:
       set id_out, file("${id_out}.fq.gz") into trimmed_fqs
       set file("*.html"), file("*.json")  into trim_log

  script:
      id_out = id.replace('.fastq.gz', '')

   """
       fastp -i $reads -o ${id_out}.fq.gz -y -l 15 -h ${id_out}.html -j ${id_out}.json
   """
}
trimmed_fqs.set { trimmed_reads_mirdeep }


////////////////////////////////////////////////
// ** - LOAD in Aedes genome, indices, and miRNA files
////////////////////////////////////////////////
// aedesgenome points to /mnt/genomes/Other/Aedes_aegypti/

geneset_stringtie = file("${aedesgenome}/annotation/geneset_h.gtf.gz")
genome_ref = file("${aedesgenome}/genome.fa")
bwa_indices = Channel.fromPath("${aedesgenome}/BWAIndex/*") //.buffer(size:8)
bowtie2_indices = Channel.fromPath("${aedesgenome}/bowtie2Index/*.bt2") //.buffer(size:8)

ae_miRNAs_mature = file(GHdata + "smRNA/miRNA/aedes_aegypti_mature_b.fasta")
ae_miRNAs_prec = file(GHdata + "smRNA/miRNA/aedes_aegypti_stemloop_b.fasta")

genome_ref.into { genome_mirdeep }


////////////////////////////////////////////////
// ** - mirDeep2 pipeline
////////////////////////////////////////////////

// Mirdeep2 mapper.pl (map to genome)
process mirDeep2_mapper_host {

    cpus large_core
    tag { id }

    input:
        set val(id), file(reads) from trimmed_reads_mirdeep
        file bowtieindex from bowtie2_indices.first()

    output:
        file("${id}_map.arf") into reads_vs_genome_arf
        set val(id), file("${id}_collapsed.fa") into reads_collapsed

    script:
        index_base = bowtie2_indices[0].toString() - ~/.\d.btw/

    """
        zcat ${reads} > ${id}.fa
        mapper.pl ${id}.fa -e -h -j -l 18 -m -p ${index_base} -s ${id}_collapsed.fa -t ${id}_map.arf -v
    """
}
reads_collapsed.into {reads_collapsed_Q; reads_collapsed_M}


// Mirdeep2 quantifier.pl (map to predefined mature/precursor seqs)
process quantifier_pl_host {

    publishDir "${output}/quantifier/${id}/", mode: 'copy'

    cpus large_core
    tag { id }

    input:
        set val(id), file(collapsed_reads) from reads_collapsed_Q

    output:
        file "*" into quantifier_out

    """
        quantifier.pl -p ${ae_miRNAs_prec} -m ${ae_miRNAs_mature} -r ${collapsed_reads} -y now
    """
}

// Mirdeep2 mirdeep2.pl
// process mirDeep2_pl {
//
//     cpus large_core
//     tag { id }
//
//     input:
//         file("genome.fa") from genome_mirdeep
//         file reads_vs_genome_arf from reads_vs_genome_arf
//         set val(id), file(collapsed_reads) from reads_collapsed_M
//
//         """
//         miRDeep2.pl ${collapsed_reads} genome.fa ${reads_vs_genome_arf} ${ae_miRNAs_mature} ${ce_miRNAs_mature} ${bm_miRNAs_prec} -P
//         """
// }
