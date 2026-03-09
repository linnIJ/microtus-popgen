#!/bin/bash
mkdir ldprune
cd ldprune
VCFFILE=../all_microtus_final.vcf.gz
TEMPID=mc_maf
FINALID=microtus_autosome_ldprune

echo "1. applying maf filter 0.01"
vcftools --gzvcf "${VCFFILE}" --maf .01 --recode --out  "${TEMPID}"
echo "2. Converting VCF to PLINK"
plink \
  --vcf "${TEMPID}".recode.vcf\
  --set-missing-var-ids @:# \
  --make-bed \
  --out "${TEMPID}" \
  --double-id \
  --vcf-half-call m \
  --chr-set 90
echo "3. Performing LD pruning on sites"
plink \
  --bfile "${TEMPID}" \
  --indep-pairwise 50 50 0.1 \
  --out "${TEMPID}"_Pruned2 \
  --chr-set 90
echo "4. Extracting LD-pruned sites"
plink \
  --bfile "${TEMPID}" \
  --extract "${TEMPID}"_Pruned2.prune.in \
  --recode vcf \
  --out "${FINALID}" \
  --chr-set 90
echo "5. Generating PLINK files from LD-Pruned Dataset"
plink \
  --vcf "${FINALID}".vcf \
  --make-bed \
  --out "${FINALID}" \
  --double-id \
  --vcf-half-call m \
  --chr-set 90

echo "6. Conducting PCA"
bcftools query -l "${FINALID}".vcf > samplenames.txt
plink --bfile "${FINALID}" --pca --out  "${FINALID}" --chr-set 90

echo "7. Conducting admixture"

for i in {2..6}; do
        admixture --cv microtus_autosome_ldprune.bed $i >log${i}.out &
done
