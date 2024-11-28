# VNCI-BreastCancer-analysis

Dependencies for  DNAdragon2.7.0
1. Python3 (tested on Python 3.6.8) with 'os', 'pandas' and 'snakemake' modules installed in SGE-cluster
2. Requires fastqc, trimmomatic, bwa, qualimap, samtools, picard, GATK versions 3.8-1-0 and >=4.1.4.1 installed

usage:

       DNAdragon2.7.0.py [-h] --conf CONF --info INFO --resource
                                RESOURCE --workDir WORKDIR --sgeLogDir
                                SGELOGDIR --maxJobPerNode MAXJOBPERNODE
                                [--keepAllFiles {T,F}] --dataType {Ex,Wg,Amp}
                                [--gvcf {T,F}] [--fixedSeqLen FIXEDSEQLEN]

one stop Exome and WGS analysis pipeline
  
optional arguments:

                -h, --help            show this help message and exit
                --conf CONF           tool config file
                --info INFO           sample sheet file
                --resource RESOURCE   refseq, 1000Genome, dbsnp data
                --workDir WORKDIR     working directory
                --sgeLogDir SGELOGDIR 
                                      sge stdout/stderr dump directory
                --maxJobPerNode MAXJOBPERNODE
                                      max job to be submitted per compute node
                --keepAllFiles {T,F}  whether to keep intermediate files (default: F)
                --dataType {Ex,Wg,Amp}
                                      Specify data type
                --gvcf {T,F}          GVCF calling (default: F)
                --fixedSeqLen FIXEDSEQLEN
                                      analysis with fixed sequence length
                Ex = Exome, Wg = WGS, Amp = Amplicon sequencing

Dependencies for somatic.sh and germline.sh
1. java-8, jdk-17
2. GATK==3.8-1-0 and >=4
3. bcftools==1.9
4. annovar==latest
5. Oncotator==1.9.9.0

Dependencies for other scripts
1. R==3.5.1 and R>=4

Script specific dependencies are listed as comments
