# Conversion of a .fam file into a regenie compatible phenotype

awk 'NR > 1 {print \$1" "\$2" "\$6}' ${FAM} > tmp_pheno.tsv
echo -e "FID IID PHENO" | cat - tmp_pheno.tsv > ${PHENO}.tsv