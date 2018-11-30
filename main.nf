#!/usr/bin/env nextflow

// Nextflow.configuration
aux=config.aux_location
data=config.data_location
output=config.output_location
GHdata=config.GHdata_location

large_core=config.large_core
small_core=config.small_core

// ** - Fetch fqs; alternative suffixes
fq_set = Channel.fromPath(data + "/*.fastq.gz")
                .map { n -> [ n.getName(), n ] }

// ** - Define paraemeters and auxillary files
adapters = file("auxillary/TruSeq3-SE.fa")
// rRNAs = file(GHdata + "smRNA/rRNA/ascaris_suum_rRNA.fasta")
// tRNAs = file(GHdata + "smRNA/tRNA/ascaris_suum_tRNA.fasta")
bm_miRNAs_mature = file(GHdata + "smRNA/miRNA/brugia_malayi_mature_b.fasta")
bm_miRNAs_prec = file(GHdata + "smRNA/miRNA/brugia_malayi_stemloop_b.fasta")
ce_miRNAs_mature = file(GHdata + "smRNA/miRNA/caenorhabditis_elegans_mature_b.fasta")
ce_miRNAs_prec = file(GHdata + "smRNA/miRNA/caenorhabditis_elegans_stemloop_b.fasta")
ae_miRNAs_mature = file(GHdata + "smRNA/miRNA/aedes_aegypti_mature_b.fasta")
ae_miRNAs_prec = file(GHdata + "smRNA/miRNA/aedes_aegypti_stemloop_b.fasta")

// ** - Fetch reference genome (fa.gz)
release="WBPS11"
species="brugia_malayi"
prjn="PRJNA10729"
prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${release}/species/${species}/${prjn}"

process fetch_parasite_ref {

    publishDir "${output}/reference/", mode: 'copy'

    output:
        file("parasite.fa.gz") into parasite_ref

    """
        echo '${prefix}'
        curl ${prefix}/${species}.${prjn}.${release}.genomic.fa.gz > parasite.fa.gz
    """
}
parasite_ref.into { parasite_bwa; parasite_bowtie; parasite_mirdeep }

// ** - Fetch host genome (fa.gz)
hosturl="https://www.vectorbase.org/download/aedes-aegypti-lvpagwgchromosomesaaegl5fagz"

// process fetch_host_ref {
//
//     publishDir "${output}/reference/", mode: 'copy'
//
//     output:
//         file("host.fa.gz") into host_ref
//
//     """
//         echo '${hosturl}'
//         wget ${hosturl} -O host.fa.gz
//     """
// }
// host_ref.into { host_bwa; host_bowtie; host_mirdeep }



//** TRIM READS
process trimmomatic {

    cpus large_core
    tag { id }
    publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*_trimout.txt'

    input:
        set val(id), file(reads) from fq_set

    output:
        set id, file(id_out) into fq_trim
        file("*_trimout.txt") into trim_log

    script:
    id_out = id.replace('.fastq.gz', '_trim.fq.gz')

    """
        trimmomatic SE -threads ${large_core} ${id} ${id_out} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 &> ${id}_trimout.txt

    """
}
fq_trim.into { fq_trim_bwa; fq_trim_bowtie; fq_trim_contam; fq_trim_mirdeepQ1; ; fq_trim_mirdeepQ2}


//INDEX GENOMES - BOWTIE
process build_bowtie_index {

    publishDir "${output}/reference/", mode: 'copy'

    cpus large_core

    input:
        file("parasite.fa.gz") from parasite_bowtie
        // file("host.fa.gz") from host_bowtie

    output:
        file "parasite_bowtie*.ebwt" into parasite_bowtie_indices
        // file "host_bowtie*.ebwt" into host_bowtie_indices

    """
        zcat parasite.fa.gz > parasite.fa
        bowtie-build parasite.fa parasite_bowtie
    """
    //        zcat host.fa.gz > host.fa
      //      bowtie-build host.fa host_bowtie
}


// Mirdeep2 quantifier.pl (map to predefined parasite mature/precursor seqs)
process quantifier_pl_parasite {

    publishDir "${output}/quantifier_parasite/", mode: 'copy'
    cpus large_core
    tag { reads }

    input:
        set val(id), file(reads) from fq_trim_mirdeepQ1

        """
        zcat ${reads} > ${fa_prefix}.fa
        quantifier.pl -p ${bm_miRNAs_prec} -m ${bm_miRNAs_mature} -r ${fa_prefix}.fa -y now
        """
}

// Mirdeep2 quantifier.pl (map to predefined host mature/precursor seqs)
process quantifier_pl_host {

    publishDir "${output}/quantifier_host/", mode: 'copy'
    cpus large_core
    tag { reads }

    input:
        set val(id), file(reads) from fq_trim_mirdeepQ2

        """
        zcat ${reads} > ${fa_prefix}.fa
        quantifier.pl -p ${ae_miRNAs_prec} -m ${ae_miRNAs_mature} -r ${fa_prefix}.fa -y now
        """
}


// // Mirdeep2 mapper.pl
// process mirDeep2_mapper {
//     cpus large_core
//     tag { id }
//
//     input:
//         set val(id), file(reads) from fq_trim_mirdeep
//         file bowtieindex from parasite_bowtie_indices.first()
//
//     output:
//         file("${fa_prefix}_parasite_map.arf") into reads_vs_parasite_genome_arf
//         file("${fa_prefix}_parasite_collapsed.fa") into reads_parasite_collapsed
//
//     script:
//         fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
//
//         """
//         zcat ${reads} > ${fa_prefix}.fa
//         mapper.pl ${fa_prefix}.fa -e -h -j -l 18 -m -p parasite_bowtie -s ${fa_prefix}_parasite_collapsed.fa -t ${fa_prefix}_parasite_map.arf -v
//         """
// }
//
//
//
// // Mirdeep2 mirdeep2.pl
// process mirDeep2_pl {
//     cpus large_core
//     tag { reads }
//
//     input:
//         file("parasite.fa.gz") from parasite_mirdeep
//         file reads_vs_parasite_genome_arf from reads_vs_parasite_genome_arf
//         file reads_parasite_collapsed from reads_parasite_collapsed
//
//         """
//         zcat parasite.fa.gz > parasite.fa
//         cat parasite.fa | awk '{print \$1}' > parasite_temp.fa
//         miRDeep2.pl ${reads_parasite_collapsed} parasite_temp.fa ${reads_vs_parasite_genome_arf} ${bm_miRNAs_mature} ${ce_miRNAs_mature} ${bm_miRNAs_prec} -P
//         """
// }


// //INDEX GENOMES - BWA
// process build_bwa_index {
//
//     publishDir "${output}/reference/", mode: 'copy'
//
//     cpus large_core
//
//     input:
//         file("parasite.fa.gz") from parasite_bwa
//         file("host.fa.gz") from host_bwa
//
//     output:
//         file "parasite.*" into bwa_parasite_indices
//         file "host.*" into bwa_host_indices
//
//     """
//         zcat parasite.fa.gz > parasite.fa
//         bwa index parasite.fa
//         zcat host.fa.gz > host.fa
//         bwa index host.fa
//     """
// }

// ALIGN TRIMMED READS TO PARASITE GENOME (BWA)
// process align {
//     publishDir "${output}/bwa_stats/", mode: 'copy'
//
//     cpus large_core
//     tag { id }
//
//     input:
//         set val(id), file(reads) from fq_trim1
//         file(parasite_bwaindex) from bwa_parasite_indices.first()
//
//     output:
//         file("bwa_parasite_align.txt") into bwa_stats
//
//     script:
//         fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
//
//         """
//         bwa aln -o 0 -n 0 -t ${large_core} parasite.fa ${reads} > ${id}.sai
//         bwa samse parasite.fa ${id}.sai ${reads} > ${id}.sam
//         samtools view -bS ${id}.sam > ${id}.unsorted.bam
//         rm *.sam
//         samtools flagstat ${id}.unsorted.bam
//         samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
//         rm *.unsorted.bam
//         samtools index -b ${id}.bam
//         samtools flagstat ${id}.bam > bwa_parasite_align.txt
//         """
// }

// // Map rRNAs and tRNAs
// process map_rRNAs_tRNAs {
//     publishDir "${output}/stats/", mode: 'copy'
//
//     cpus large_core
//     tag { reads }
//
//     input:
//         file reads from fq_trim2
//         file rRNA_fa from rRNAs
//         file tRNA_fa from tRNAs
//
//     output:
//         file("bwa_rRNA_align.txt") into bwa_rRNA_alignstats
//         file("bwa_tRNA_align.txt") into bwa_tRNA_alignstats
//
//     script:
//         fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
//
//         """
//         bwa index ${rRNA_fa}
//
//         bwa aln -o 0 -n 0 -t ${large_core} ${rRNA_fa} ${reads} > ${fa_prefix}.sai
//         bwa samse ${rRNA_fa} ${fa_prefix}.sai ${reads} > ${fa_prefix}.sam
//         samtools view -bS ${fa_prefix}.sam > ${fa_prefix}.unsorted.bam
//         samtools flagstat ${fa_prefix}.unsorted.bam
//         samtools sort -@ ${large_core} -o ${fa_prefix}_rRNA.bam ${fa_prefix}.unsorted.bam
//         samtools index -b ${fa_prefix}_rRNA.bam
//         samtools flagstat ${fa_prefix}_rRNA.bam > bwa_rRNA_align.txt
//
//         bwa index ${tRNA_fa}
//
//         bwa aln -o 0 -n 0 -t ${large_core} ${tRNA_fa} ${reads} > ${fa_prefix}.sai
//         bwa samse ${tRNA_fa} ${fa_prefix}.sai ${reads} > ${fa_prefix}.sam
//         samtools view -bS ${fa_prefix}.sam > ${fa_prefix}.unsorted.bam
//         samtools flagstat ${fa_prefix}.unsorted.bam
//         samtools sort -@ ${large_core} -o ${fa_prefix}_tRNA.bam ${fa_prefix}.unsorted.bam
//         samtools index -b ${fa_prefix}_tRNA.bam
//         samtools flagstat ${fa_prefix}_tRNA.bam > bwa_tRNA_align.txt
//
//         """
// }
