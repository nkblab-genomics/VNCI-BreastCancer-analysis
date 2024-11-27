
#!/bin/bash
Thisdir=</Directory for analysis>
cd "$Thisdir"

# Dependency: java-8, gatk-4.0.11.0
# For this project chrY is excluded as all individuals are female

REF=</Reference fasta file>
GenomeADb37=af-only-gnomad.raw.sites.b37.vcf.gz

Sample_id="$1"
Sample_Normal="$2"
BAM_Normal="$3"
Sample_Tumor="$4"
BAM_Tumor="$5"

mkdir -p TMP_DIR

gatk --java-options "-Xms14g -Xmx16g" Mutect2 \
-R $REF \
-I $BAM_Tumor \
-tumor $Sample_Tumor \
-I $BAM_Normal \
-normal $Sample_Normal \
--germline-resource $GenomeADb37 \
-O "$Sample_id".mutect2.vcf.gz \
--annotation StrandBiasBySample \
--native-pair-hmm-threads 4 \
--tmp-dir TMP_DIR \
--panel-of-normals </panel of normal vcf.gz> \
-XL chrY

$GATKHOME/gatk --java-options "-Xms14g -Xmx16g" FilterMutectCalls \
-R "$REF" \
-V "$Sample_id".mutect2.vcf.gz \
-O "$Sample_id".mutect2_filter.vcf.gz \
-XL chrY
