#!/bin/bash
set -euo pipefail

echo "1: remove duplicate and peromyscus, then keep  variant sites"
echo MVZ_Mamm_227885_Mca69 > dup.txt
echo MVZ_Mamm_207423_Mca9 > peromyscus.txt
cat dup.txt peromyscus.txt > dup_peromyscus.txt
bcftools view ../41-Microtus_raw.vcf.gz --samples-file ^dup_peromyscus.txt -Ou | \
bcftools +fill-tags -- -t AN,AC,AF | \
bcftools view -i 'MAF>0' -Oz -o all_microtus.vcf.gz
rm dup.txt
rm peromyscus.txt
rm dup_peromyscus.txt

echo "2: set sites with QUAL/GQ less than 20 to missing"
bcftools filter -S . -e 'QUAL<20' all_microtus.vcf.gz -Ou | \
bcftools filter -S . -e 'FORMAT/GQ<20' -Oz -o microtus_qualityfiltered.vcf.gz
bcftools index -t microtus_qualityfiltered.vcf.gz

echo "3: Annotate vcf with numerical chromosomes, then recalculate allele frequencies"
bcftools annotate --rename-chrs chrommap.txt microtus_qualityfiltered.vcf.gz -Oz -o temp_annotated.vcf.gz
bcftools index -t temp_annotated.vcf.gz
bcftools +fill-tags temp_annotated.vcf.gz -Oz -o microtus_annotated.vcf.gz -- -t AN,AC,AF
bcftools index -t microtus_annotated.vcf.gz

echo "4: Subset to X and autosomes, excluding variants missing in 25% samples, mean depth lower than 5 or greater than 12, and multiallelic sites"
vcftools --gzvcf microtus_annotated.vcf.gz --max-missing 0.75  --max-alleles 2\
  --chr 1 --chr 2 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 \
  --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 \
  --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27  --min-meanDP 5 --max-meanDP 15\
  --recode --out all_microtus_cleaned

bgzip -c all_microtus_cleaned.recode.vcf> all_microtus_cleaned.recode.vcf.gz
tabix -p vcf all_microtus_cleaned.recode.vcf.gz

echo "5: subset to bedfile callable regions"
bcftools view -R 41-Microtus_callable_renamed.bed all_microtus_cleaned.recode.vcf.gz -o all_microtus_final.vcf.gz
bcftools index -t all_microtus_final.vcf.gz






