
#!/bin/bash
Thisdir=</Directory for analysis>
cd "$Thisdir"

# Dependency: java-8, gatk-3.8-1-0 and >=4, bcftools-1.9, CombineSomatic-0.2, oncotator==1.9.9.0
# For this project chrY is excluded as all individuals are female

REF=</Reference fasta file>
GenomeADb37=af-only-gnomad.raw.sites.b37.vcf.gz

Sample_id="$1"
Sample_Normal="$2"
BAM_Normal="$3"
Sample_Tumor="$4"
BAM_Tumor="$5"

mkdir -p TMP_DIR

# run using GATK>=4
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

gatk --java-options "-Xms14g -Xmx16g" FilterMutectCalls \
-R "$REF" \
-V "$Sample_id".mutect2.vcf.gz \
-O "$Sample_id".mutect2_filter.vcf.gz \
-XL chrY

# run using GATK==3.8-1-0
java -jar -Xmx4g GenomeAnalysisTK-3.8-1-0/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R $REF -V "$Sample_id".mutect2_filter.vcf.gz -o "$Thisdir"/"$Sample_id"_temp1.vcf
bcftools norm --check-ref e --fasta-ref $REF --multiallelics '-both' --output-type v "$Sample_id"_temp1.vcf > "$Thisdir"/"$Sample_id"_mutect2_normalised.vcf
rm -f "$Thisdir"/"$Sample_id"_temp1.vcf*

# fetch depth information from vcf
python CombineSomatic-0.2.py --config config.CombineSomatic.json --sample "$Sample_id" --variantCalls Mutect2:"$Sample_id"_mutect2_normalised.vcf --out "$Sample_id".depth.mutect2.tsv

# annotation through oncotator==1.9.9.0
grep "^""#" -v "$Thisdir"/"$Sample_id"_mutect2_normalised.vcf | cut -f1-5 | awk '{if(length($4)>1) {print $1"\t"$2+1"\t"$2+length($4)-1"\t"substr($4,2,length($4))"\t""-""\t"$0} else if(length($5)>1) {print $1"\t"$2"\t"$2"\t""-""\t"substr($5,2,length($5))"\t"$0} else {print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$0}}' > "$Sample_id"_m2_v0.txt
echo -e "chr\tstart\tend\tref_allele\talt_allele" > "$Sample_id"_m2.txt
cut -f1-5 "$Sample_id"_m2_v0.txt >> "$Sample_id"_m2.txt
Oncotator -v --db-dir </path for oncotator_v1_ds_April052016> -i MAFLITE -o TCGAMAF --skip-no-alt "$Sample_id"_m2.txt "$Sample_id"_m2.somatic.oncotator.maf hg19

# "$Sample_id".depth.mutect2.tsv and "$Sample_id"_m2.somatic.oncotator.maf can be combined based on chr, start, end, ref_allele and alt_allele for final filtration described in Ghosh et al., manuscript.
