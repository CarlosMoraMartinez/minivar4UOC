#Create a new conda environment with python 3.10
conda create -n  minivar4UOC_env python=3.10 numpy pandas 

#Activate the environment
conda activate minivar4UOC_env

#Some bioinformatics basic software to work with bam, vcf, bed and fastq files
conda install -c bioconda bedtools samtools seqtk bcftools vt bedops
conda install -c bioconda tabix

#Trim Galore for trimming reads. Installs FastQC and cutadapt
conda install -c bioconda trim-galore

#Also, let's try MultiQC for aggregated quality control and fastp
conda install -c bioconda fastp # not
conda install -c bioconda multiqc

#bwa for alignment
conda install -c bioconda bwa

#Some variant callers:
conda install -c bioconda freebayes
conda install -c bioconda deepvariant ##failed
conda install -c bioconda vardict

 #Deeptools for plots and downstream analysis
conda install -c bioconda deeptools #failed

#Filter and annotate variants
#conda install -c bioconda snpeff
#conda install -c "bioconda/label/cf201901" ensembl-vep

#Download snpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
java -jar snpEff.jar download GRCh37.75