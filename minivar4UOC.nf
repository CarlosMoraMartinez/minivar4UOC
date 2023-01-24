#!/usr/bin/env nextflow


process trimReads{
  cpus params.resources.trim.cpus
  memory params.resources.trim.mem
  publishDir "$results_dir/01_trimmed", mode: 'symlink'
  input:
  path fastq
  val quality
  val min_length
  
  output:
  tuple(path("*trimmed.fq.gz"), path('*_trimming_report.txt'))

  shell:
  '''
  trim_galore -q !{quality} --length !{min_length} --gzip !{fastq}
  '''
}

process getFastQC{
  cpus params.resources.fastqc.cpus
  memory params.resources.fastqc.mem
  publishDir "$results_dir/02_fastqc", mode: 'copy'
  input:
  path fastq
  
  output:
  path '*.html'

  shell:
  '''
  fastqc -q !{fastq}  
  '''
}

process buildIndex {
  cpus params.resources.index.cpus
  memory params.resources.index.mem
  publishDir "$results_dir/00_index", mode: 'copy'

  input:
  path fasta_file
  
  output:
  tuple path("${fasta_file}"), path("${fasta_file}*")


  shell:
  '''
    bwa index !{fasta_file}
  '''
}

process buildFaiIndex {
  cpus params.resources.index.cpus
  memory params.resources.index.mem
  publishDir "$results_dir/00_index", mode: 'copy'

  input:
  path fasta_file
  
  output:
  tuple path("*.fa"), path("*.fai")


  shell:
  '''   
    fastagenome=$(basename -s .gz !{fasta_file})
    zcat !{fasta_file} > $fastagenome
    samtools faidx $fastagenome
    
  '''
}

process alignReads {
  cpus params.resources.alignment.cpus
  memory params.resources.alignment.mem
  publishDir "$results_dir/03_align", mode: 'symlink'

  input:
    val genome
    val illumina_reads 

  output:
    tuple(path('*_raw.bam'),path('*_raw.bam.bai'))

  shell:
    '''
    outfile=$(basename -s .fq.gz !{illumina_reads})'_raw.bam'

    bwa mem -t !{params.resources.alignment.cpus} !{params.alignment.bwaparams} !{genome}  !{illumina_reads} | samtools sort -o $outfile
    samtools index $outfile
    '''
}


process filterAlignment{
  cpus params.resources.samtoolsfilter.cpus
  memory params.resources.samtoolsfilter.mem
  publishDir "$results_dir/04_filter-F$exclude_flag_F-q$mapq", mode: 'symlink'

  input:
    path bam
    val mapq
    val include_flag_f
    val exclude_flag_F

  output:
    tuple(path('*.filt.bam'),path('*.filt.bam.bai'))

  shell:
  '''
  outfile=$(basename -s _raw.bam !{bam})
  filtbam=$outfile'.filt.bam'

  samtools view -@ !{params.resources.samtoolsfilter.cpus} -b -F !{exclude_flag_F} -f !{include_flag_f} -q !{mapq} !{bam} -o $filtbam
  samtools index $filtbam
  '''
}

process removeDups{
  cpus params.resources.samtoolsfilter.cpus
  memory params.resources.samtoolsfilter.mem
  publishDir "$results_dir/05_rmdup", mode: 'symlink'

  input:
    path bam

  output:
    tuple(path('*.rmdup.bam'), path('*.rmdup.bam.bai'))

  shell:
  '''
  outfile=$(basename -s .bam !{bam})
  filtbam=$outfile'.rmdup.bam'

  samtools rmdup -s !{bam} $filtbam
  samtools index $filtbam
  '''
}

process callVariants {
  cpus params.resources.alignment.cpus
  memory params.resources.alignment.mem
  publishDir "$results_dir/06_vars", mode: 'symlink'

  input:
    val genome
    val bam 

  output:
    tuple(path('*.vcf.gz'),path('*.vcf.gz.tbi'))

  shell:
    '''
    outfile=$(basename -s .bam !{bam})'.vcf'
    freebayes -f !{genome} !{bam} > $outfile
    bgzip $outfile
    tabix $outfile'.gz'
    '''
}

process annotateSnpEff{
  cpus params.resources.cov.cpus
  memory params.resources.cov.mem
  publishDir "$results_dir/07_annot", mode: 'symlink'

  input:
    val vcf
    val execdir
    val snpeff_db
    val external_db


  output:
    tuple(path("*.vcf.gz"),path('*.vcf.gz.tbi'), path('*.html'), path('*.txt'))

  shell:
    '''
    basename=$(basename -s '.vcf.gz' !{vcf})
    outname=$basename'.ann.vcf'
    java -Xmx8g -jar !{execdir}/snpEff.jar -v -stats $basename'.html' !{snpeff_db} !{vcf} | \
    java -Xmx8g -jar !{execdir}/SnpSift.jar annotate -v !{external_db} > $outname
    bgzip $outname
    tabix $outname'.gz'   
    '''
}

process filterVars{
  cpus params.resources.cov.cpus
  memory params.resources.cov.mem
  publishDir "$results_dir/08_filtvcf_dp$params.filter.mindp-freq$params.filtvcf.filter_freq-ann$params.filtvcf.filter_ann", mode: 'symlink'

  input:
    tuple(val(vcf), val(tbi))
    val mindp
    val filter_freq
    val af_exac
    val af_esp
    val af_tgp
    val filter_ann


  output:
    tuple(path("*filt.vcf.gz"),path('*filt.vcf.gz.tbi'))

  shell:
    '''
    basename=$(basename -s '.vcf.gz' !{vcf})
    outname=$basename'.filt.vcf'
    filter_annot=!{filter_ann}

    if [ "$filter_annot" = "" ]; then
      typefilt=""
    else
      typefilt=' && INFO/ANN[*]~"'"$filter_annot"'"'
    fi

    if [ "!{filter_freq}" = "false" ]; then
      freq_filter=''
    else
      freq_filter='(INFO/AF_EXAC < !{af_exac} || INFO/AF_ESP < !{af_esp} || INFO/AF_TGP < !{af_tgp}) && '
    fi

    bcftools filter -i "$freq_filter"'INFO/DP > !{mindp}'"$typefilt" !{vcf} > $outname 

    bgzip $outname
    tabix $outname'.gz'   
    '''
}


process calcCoverage{
  cpus params.resources.cov.cpus
  memory params.resources.cov.mem
  publishDir "$results_dir/09_coverage", mode: 'symlink'

  input:
    val bam
    val genome_sizes
    val bed

  output:
    path("*.tsv")

  shell:
    '''
    outname=$(basename -s '.bam' !{bam})'.tsv'
    bedtools coverage -a !{bed} -b !{bam} -sorted -mean -g !{genome_sizes} > $outname
    '''
}


process reportExonGaps{
  cpus params.resources.cov.cpus
  memory params.resources.cov.mem
  publishDir "$results_dir/10_gaps", mode: 'symlink'

  input:
    val tsv
    val dp

  output:
    path("*.gaps.bed")

  shell:
    '''
    outname=$(basename -s '.tsv' !{tsv})'.gaps.bed'
    dplim=!{dp}
    awk -F '\t' -v dp=$dplim '{ if ($7 < dp) { print } }' !{tsv} > $outname
    '''
}

workflow {

  //Build genome index if needed
  if(params.index.do){
    buildIndex(params.genomefasta)
    ch_index = buildIndex.out
        .view{"Genome Index: $it"}
        .map{it -> it[0]}
  }else{
    ch_index = Channel.fromPath(params.genomeindex)
  }

  //Build .fai index for variant calling if needed
  if(params.indexFai.do){
    buildFaiIndex(params.genomefasta)
    ch_indexfai = buildFaiIndex.out
        .view{"Genome Index Fai: $it"}
        .map{it -> it[0]}
  }else{
    ch_indexfai = Channel.fromPath(params.genomeindexfai)
  }


  //Trim reads 
   ch_reads = Channel.fromPath(params.reads)
   if(params.fastqc.do){
    ch_fastqc = ch_reads 
   }
   trimReads(ch_reads, params.trim.quality, params.trim.min_length)
   ch_reads = trimReads.out
    .view{"Trimmed fastq: $it"}
    .map{it -> it[0]}

  //Calculate FastQC
   if(params.fastqc.do){
    ch_fastqc = ch_fastqc.concat(ch_reads)
    .view{"FastQC input: $it"}
    getFastQC(ch_fastqc)
    .view{"FastQC output: $it"}
   }

  //Align reads
   alignReads(ch_index, ch_reads)
   ch_bam = alignReads.out
   .map{it -> it[0]}
   .view{"Aligned reads raw: $it"}
   ch_cov = ch_bam

  //Filter by mapping quality and flags
  filterAlignment(
    ch_bam,
    params.filterAlignment.mapq,
    params.filterAlignment.include_flag_f,
    params.filterAlignment.exclude_flag_F
  )
 
  ch_filtered = filterAlignment.out
  .map{it -> it[0]}
  .view{ "Filtered reads: $it" }
  ch_cov = ch_cov.concat(ch_filtered)

  //Remove reads
  removeDups(ch_filtered)
  ch_filtered = removeDups.out
  .map{it -> it[0]}
  .view{ "Rmdup result reads: $it" }
  ch_cov = ch_cov.concat(ch_filtered)

  //Call variants
  callVariants(ch_indexfai, ch_filtered)
  ch_variants = callVariants.out
  .map{it -> it[0]}
  .view{ "Variants: $it" }

  //Annotate vcf
  annotateSnpEff(ch_variants, params.snpeff.execdir, params.snpeff.db, params.snpeff.external_db)
  ch_variants = annotateSnpEff.out 
   .view{ "Annotated Variants: $it" }
   .map{it -> [it[0], it[1]] }

  //Filter variants
  filterVars(ch_variants,
    params.filter.mindp,
    params.filtvcf.filter_freq,
    params.filtvcf.af_exac,
    params.filtvcf.af_esp,
    params.filtvcf.af_tgp,
    params.filtvcf.filter_ann
)
  ch_variants = filterVars.out 
   .view{ "Filtered Variants: $it" }
    
  //Calculated coverage for CDS regions
  calcCoverage(ch_cov, params.genomesizes, params.coverage.bed)
  ch_cov = calcCoverage.out
  .view{ "Coverage: $it" }

  //Report regions in which variants will not be detected
  reportExonGaps(ch_cov, params.filter.mindp)
  ch_gaps = reportExonGaps.out
    .view{ "Gaps: $it" }

}