#!/bin/bash
# Dependency: java-8, GATK-4, bcftools-1.9, tabix, bgzip, annovar

Thisdir=</Path for analysis>
cd "$Thisdir"

REF=</Reference genome fasta>
GATK_JAR=gatk-package-4.x.x.x-local.jar 
humandb_path=</annovar humandb path hg19>

# add germline g.vcf.gz for all samples, the example is provided for two samples.
java -Xms22g -Xmx23g -jar "$GATK_JAR" CombineGVCFs \
 -R "$REF" \
 --variant </sample1.g.vcf.gz> \
 --variant </sample2.g.vcf.gz> \
 -O cohort.joint.g.vcf.gz

java -Xms22g -Xmx23g -jar "$GATK_JAR" GenotypeGVCFs \
 -R "$REF" \
 -V cohort.joint.g.vcf.gz \
 -O cohort.joint.vcf.gz

# annotation
bcftools view -M 2 -Oz -o cohort.joint.bi_allele.vcf.gz cohort.joint.vcf.gz
tabix -p vcf cohort.joint.bi_allele.vcf.gz

gatk VariantAnnotator -R "$REF" \
 --variant cohort.joint.bi_allele.vcf.gz \
 --resource:1000G 1000G_phase3_v4_20130502.sites.chr.vcf.gz \
 --resource:gnomad af-only-gnomad.raw.sites.b37.vcf.gz \
 --resource:clinvar clinvar_20241027.chr.vcf.gz \
 -E 1000G.AF  -E gnomad.AF  -E clinvar.CLNSIG  -E clinvar.GENEINFO  -E clinvar.MC \
 --resource-allele-concordance \
 -O cohort.joint.bi_allele.1000G_genomad_clinvar.vcf.gz \
 --disable-sequence-dictionary-validation true \
 --java-options "-Xms12g -Xmx64g"

bgzip -c cohort.joint.bi_allele.1000G_genomad_clinvar.vcf.gz > cohort.joint.bi_allele.1000G_genomad_clinvar.vcf

annovar/table_annovar.pl cohort.joint.bi_allele.1000G_genomad_clinvar.vcf \
 "$humandb_path" -buildver hg19 \
 -out cohort.joint.bi_allele.1000G_genomad_clinvar.annovar \
 -remove -protocol refGene,cytoBand,dbnsfp30a,intervar_20170202,avsnp147 \
 -operation gx,r,f,f,f -nastring . -polish -vcfinput
