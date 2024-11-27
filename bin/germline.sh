#!/bin/bash
# Dependency: java-8, GATK-4

Thisdir=</Path for analysis>
cd "$Thisdir"

REF=</Reference genome fasta>
GATK_JAR=gatk-package-4.x.x.x-local.jar 

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
