# SAMPLE, GEN, GENTIC_MAP, HAPS, LEGEND, STRAND_INFO
# imputed.gen, imputed.gen_warnings

chrXflags=''
if [[ ${CHR} == *X* ]]
  then
    chrXflags="-sample_g ${SAMPLE} -chrX"
    if [[ ${CHR} != *NONPAR* ]]
    then
        chrXflags="\$chrXflags -Xpar"
    fi
fi

grep "^${CHR} " ${STRAND_INFO} | cut -d' ' -f2,3 >strand_g_file
grep "^${CHR} " ${GEN} | awk '\$3 >= ${START}'| awk '\$3 <= ${END}' >region.gen

# impute only on the screened SNPs
## only genotyped snps are included in the panel
## use 503 samples for the reference with european origin
impute2 \
-g region.gen \
-m ${REFERENCE}/genetic_map_chr${CHR}_combined_b37.txt \
-int ${START} ${END} \
-h ${REFERENCE}/1000GP_Phase3_chr${CHR}.hap.gz \
-l ${REFERENCE}/1000GP_Phase3_chr${CHR}.legend.gz \
-Ne 20000 \
-verbose \
\$chrXflags \
-k_hap 503 \
-os 2 -o imputed.gen \
-strand_g strand_g_file
