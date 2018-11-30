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
fq_trim.into { fq_trim_bwa; fq_trim_bowtie; fq_trim_contam; fq_trim_mirdeep }


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

// Mirdeep2 mapper.pl
process mirDeep2_mapper {
    cpus large_core
    tag { id }

    input:
        set val(id), file(reads) from fq_trim_mirdeep
        file bowtieindex from parasite_bowtie_indices.first()

    output:
        file("${fa_prefix}_parasite_map.arf") into reads_vs_parasite_genome_arf
        file("${fa_prefix}_parasite_collapsed.fa") into reads_parasite_collapsed

    script:
        fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/

        """
        zcat ${reads} > ${fa_prefix}.fa
        mapper.pl ${fa_prefix}.fa -e -h -j -l 18 -m -p parasite_bowtie -s ${fa_prefix}_parasite_collapsed.fa -t ${fa_prefix}_parasite_map.arf -v
        """
}


// Mirdeep2 mirdeep2.pl
process mirDeep2_pl {
    cpus large_core
    tag { reads }

    input:
        file("parasite.fa.gz") from parasite_mirdeep
        file reads_vs_parasite_genome_arf from reads_vs_parasite_genome_arf
        file reads_parasite_collapsed from reads_parasite_collapsed

        """
        zcat parasite.fa.gz > parasite.fa
        cat parasite.fa | awk '{print \$1}' > parasite_temp.fa
        miRDeep2.pl ${reads_parasite_collapsed} parasite_temp.fa ${reads_vs_parasite_genome_arf} ${bm_miRNAs_mature} ${ce_miRNAs_mature} ${bm_miRNAs_prec} -P
        """
}


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




//** - ALIGNMENT AND STRINGTIE (combined)
// process align_stringtie {

//     publishDir "${output}/expression", mode: 'copy'

//     cpus large_core

//     tag { id }

//     input:
//         set val(id), file(forward), file(reverse) from read_pairs
//         file("geneset.gtf.gz") from mm_gtf
//         //file("genome_tran.5.ht2"), file("genome_tran.3.ht2"), file("genome_tran.4.ht2"), file("genome_tran.6.ht2"), file("genome_tran.1.ht2"), file("genome_tran.8.ht2"), file("genome_tran.2.ht2"), file("genome_tran.7.ht2") from hs2_indices

//     output:
//         file "${id}.hisat2_log.txt" into alignment_logs
//         file("${id}/*") into stringtie_exp

//     //script:
//      //   index_base = hs2_indices[0].toString() - ~/.\d.ht2/

//     """
//         hisat2 -p ${large_core} -x '/home/BIOTECH/zamanian/GitHub/AsELV_RNAseq-nf/data/reference/grcm38_tran/genome_tran' -1 ${forward} -2 ${reverse} -S ${id}.sam --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" 2> ${id}.hisat2_log.txt
//         samtools view -bS ${id}.sam > ${id}.unsorted.bam
//         rm *.sam
//         samtools flagstat ${id}.unsorted.bam
//         samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
//         rm *.unsorted.bam
//         samtools index -b ${id}.bam
//         zcat geneset.gtf.gz > geneset.gtf
//         stringtie ${id}.bam -p ${large_core} -G geneset.gtf -A ${id}/${id}_abund.tab -e -B -o ${id}/${id}_expressed.gtf
//         rm *.bam
//         rm *.bam.bai
//         rm *.gtf
//     """
// }



//comment out until all else finished
// prepDE = file("${aux}/scripts/prepDE.py")

// process stringtie_table_counts {

//     echo true

//     publishDir "${output}/diffexp", mode: 'copy'

//     cpus small_core

//     output:
//         file ("gene_count_matrix.csv") into gene_count_matrix
//         file ("transcript_count_matrix.csv") into transcript_count_matrix

//     """
//         python ${prepDE} -i ${output}/expression -l 140 -g gene_count_matrix.csv -t transcript_count_matrix.csv

//     """
// }








// // // // // // // // // // // // // // // // // IGNORE BELOW

// ** - Recurse through subdirectories to get all fastqs
// fq_set = Channel.fromPath(data + "fq/*.fastq.gz")
//                 .map { n -> [ n.getName(), n ] }

// SKIP TRIMMING (READS ARE ALREADY TRIMMED)
// process trim {

//     tag { fq_id }

//     publishDir "${data}/fq_trim/", mode: 'move'

//     input:
//         set fq_id, file(forward), file(reverse) from read_pairs

//     output:
//         set file("${fq_id}_1P.fq.gz"), file("${fq_id}_2P.fq.gz") into trim_output
//         //file "${fq_id}.trim_log.txt" into trim_logs

//     """
//     trimmomatic PE -threads ${large_core} $forward $reverse -baseout ${fq_id}.fq.gz ILLUMINACLIP:/home/linuxbrew/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:80:10 MINLEN:75
//     rm ${fq_id}_1U.fq.gz
//     rm ${fq_id}_2U.fq.gz
//     """

// }


// process trimmomatic {

//     cpus small_core

//     tag { name }

//     input:
//         set val(name), file(reads) from fq_set

//     output:
//         file(name_out) into trimmed_reads

//     script:
//     name_out = name.replace('.fastq.gz', '_trim.fq.gz')

//     """
//         trimmomatic SE -phred33 -threads ${small_core} ${reads} ${name_out} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
//     """
// }







// process stringtie_counts {

//     publishDir "output/expression", mode: 'copy'

//     cpus small_core

//     tag { srid }

//     input:
//         set val(srid), file(bam), file(bai) from hisat2_bams
//         file("geneset.gtf.gz") from geneset_stringtie.first()

//     output:
//         file("${srid}/*") into stringtie_exp

//     """
//         zcat geneset.gtf.gz > geneset.gtf
//         stringtie -p ${small_core} -G geneset.gtf -A ${srid}/${srid}_abund.tab -e -B -o ${srid}/${srid}_expressed.gtf ${bam}
//     """
// }



// // ** - ALIGNMENT
// process align {

//     cpus small_core

//     tag { srid }

//     input:
//         set val(srid), file(forward), file(reverse) from read_pairs
//         file hs2_indices from hs2_indices.first()

//     output:
//         set val(srid), file("${srid}.bam"), file("${srid}.bam.bai") into hisat2_bams
//         file "${srid}.hisat2_log.txt" into alignment_logs

//     script:
//         index_base = hs2_indices[0].toString() - ~/.\d.ht2/

//     """
//         hisat2 -p ${small_core} -x $index_base -1 ${forward} -2 ${reverse} -S ${srid}.sam --rg-id "${srid}" --rg "SM:${srid}" --rg "PL:ILLUMINA" 2> ${srid}.hisat2_log.txt
//         samtools view -bS ${srid}.sam > ${srid}.unsorted.bam
//         samtools flagstat ${srid}.unsorted.bam
//         samtools sort -@ ${small_core} -o ${srid}.bam ${srid}.unsorted.bam
//         samtools index -b ${srid}.bam
//         rm *sam
//         rm *unsorted.bam

//     """
// }



// process stringtie_counts {

//     publishDir "output/expression", mode: 'copy'

//     cpus small_core

//     tag { srid }

//     input:
//         set val(srid), file(bam), file(bai) from hisat2_bams
//         file("geneset.gtf.gz") from geneset_stringtie.first()

//     output:
//         file("${srid}/*") into stringtie_exp

//     """
//         zcat geneset.gtf.gz > geneset.gtf
//         stringtie -p ${small_core} -G geneset.gtf -A ${srid}/${srid}_abund.tab -e -B -o ${srid}/${srid}_expressed.gtf ${bam}
//     """
// }



// prepDE = file("auxillary/scripts/prepDE.py")

// process stringtie_table_counts {

//     echo true

//     publishDir "output/diffexp", mode: 'copy'

//     cpus small_core

//     tag { sample_id }

//     input:
//         val(sample_file) from stringtie_exp.toSortedList()

//     output:
//         file ("gene_count_matrix.csv") into gene_count_matrix
//         file ("transcript_count_matrix.csv") into transcript_count_matrix

//     """
//         for i in ${sample_file.flatten().join(" ")}; do
//             bn=`basename \${i}`
//             full_path=`dirname \${i}`
//             sample_name=\${full_path##*/}
//             echo "\${sample_name} \${i}"
//             mkdir -p expression/\${sample_name}
//             ln -s \${i} expression/\${sample_name}/\${bn}
//         done;
//         python ${prepDE} -i expression -l 50 -g gene_count_matrix.csv -t transcript_count_matrix.csv

//     """
// }
