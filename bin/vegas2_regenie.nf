#!/usr/bin/env nextflow

// Localisation of the output
params.out = '.'
// Bed file containing genotype
params.bed = ''
// Version of gencode to use for retrieving the gene id list
params.gencode = 28
params.genome = '38'
// Additional parameters for VEGAS2
params.vegas_params = ''
// If provided, regenie isn't used to compute snp level association scores
params.snp_association = ''
// If provided, used by VEGAS2 as a control for LD when computing gene level association scores
params.ld_controls = ''
// Buffer for gff2bed
params.buffer = 0
// Phenotype file to be used by regenie
params.phenotype = ''
// Covar file to be used by regenie
params.covar = 'NO_COV_FILE'
// regenie genotype block size for step 1
params.bsize_s1 = 1000
// regenie genotype block size for step 2
params.bsize_s2 = 400
// Additional alphanumerical parameters for regenie step 1
params.regenie_params_s1 = ''
// Additional alphanumerical parameters for regenie step 2
params.regenie_params_s2 = ''
// Mapping of the main regenie's parameters linking a file
params.extract_snp = 'NO_EXTRACT_FILE'
params.exclude_snp = 'NO_EXCLUDE_FILE'
params.keep_indiv = 'NO_KEEP_FILE'
params.remove_indiv = 'NO_REMOVE_FILE'
// output folder
params.output = '.'
// subfolder for regenie indermediate files
params.regenie_folder = 'regenie_run'

// Nextflow "consume" params variables once it's used, so retrieving
// csorresponding file for use later
regenie_folder = "${params.regenie_folder}"
output_folder = file("${params.output}")
// Plink works with a incomplete path for bed, bim, fam but files have to be 
// specified as input in nextflow pipelines
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

// Multiple process use this files so dedicating a channel
// Note: not doing covar here as it is used only in one process
extract_path = file("${params.extract_snp}")
exclude_path = file("${params.exclude_snp}")
keep_path = file("${params.keep_indiv}")
remove_path = file("${params.remove_indiv}")

// Checking double dashes have been escaped for regenie_params_sX
if (params.regenie_params_s1 == true) {
    exit 1, "Escape of `--` for regenie_params_s1 content required. Ex: `\\-\\-cc12`"
}
if (params.regenie_params_s2 == true) {
    exit 1, "Escape of `--` for regenie_params_s2 content required. Ex: `\\-\\-cc12`"
}

//////////////////////////////////////////////
///            CREATE GENE LIST            ///
//////////////////////////////////////////////
process download_gencode {

    input:
        val GENCODE_VERSION from params.gencode
        val GRCH_VERSION from params.genome

    output:
        file 'gff3' into gff
    
    script:
    template 'dbs/gencode.sh'

}

process make_glist {

    input:
        file gff
        val BUFF from params.buffer

    output:
        file 'glist_ensembl' into glist_vegas, glist_chrs

    """
    awk '\$3 == "gene"' $gff >genes.gff
    gff2bed < genes.gff | cut -f1-4 | sed 's/\\.[^\\t]\\+\$//' | sed 's/^chr//' >tmp
    awk '{\$2 = \$2 - ${BUFF}; \$3 = \$3 + ${BUFF}} 1' tmp | awk '\$2 < 0 {\$2 = 0} 1' >buffered_genes
    sed 's/^XY/25/' buffered_genes | sed 's/^X/23/' | sed 's/^Y/24/' | sed 's/^M/26/' | awk '\$1 <= 24' >glist_ensembl
    """

}

chromosomes = glist_chrs
    .splitCsv(header: false, sep: ' ')
    .map { row -> row[0] }
    .unique()

//////////////////////////////////////////////
///          SNP ASSOCIATION TEST          ///
//////////////////////////////////////////////

extract = extract_path.name != 'NO_EXTRACT_FILE' ? "--extract ${extract_path}" : ''
exclude = exclude_path.name != 'NO_EXCLUDE_FILE' ? "--exclude ${exclude_path}" : ''
keep = keep_path.name != 'NO_KEEP_FILE' ? "--keep ${keep_path}" : ''
remove = remove_path.name != 'NO_REMOVE_FILE' ? "--remove ${remove_path}" : ''

// If no snp association is provided
if (params.snp_association == ''){

    // Getting the covariable file into a chanel
    covar_path = file("${params.covar}")

    process compute_association_score {

        publishDir "${output_folder}/${regenie_folder}", overwrite: true, mode: "copy"

        input:
            file BED from bed
            file bim
            file fam
            file PHENO from file("${params.phenotype}")
            file OUTPUT from output_folder
            val  BSIZE_S1 from params.bsize_s1
            val  BSIZE_S2 from params.bsize_s2
            val  REGENIE_PARAMS_S1 from params.regenie_params_s1
            val  REGENIE_PARAMS_S2 from params.regenie_params_s2
            path covar_path
            path extract_path
            path exclude_path
            path keep_path
            path remove_path
            
        output:
            file 'regenie_output_s2_*.regenie' into regenie_output_ch
            file 'regenie_output_s*'

        script:
        covar = covar_path.name != 'NO_COV_FILE' ? "--covarFile ${covar_path}" : ''
        """
        regenie \\
            --step 1 \\
            --bed ${BED.baseName} \\
            --phenoFile ${PHENO} \\
            $covar \\
            $extract \\
            $exclude \\
            $keep \\
            $remove \\
            --bt \\
            --bsize ${BSIZE_S1} \\
            --lowmem \\
            --lowmem-prefix tmp_lowmem_regenie_s1 \\
            --out regenie_output_s1 \\
            ${REGENIE_PARAMS_S1}

        regenie \\
            --step 2 \\
            --bed ${BED.baseName} \\
            --phenoFile ${PHENO} \\
            $covar \\
            $extract \\
            $exclude \\
            $keep \\
            $remove \\
            --bt \\
            --bsize ${BSIZE_S2} \\
            --pred regenie_output_s1_pred.list \\
            --print-pheno \\
            --out regenie_output_s2 \\
            ${REGENIE_PARAMS_S2}
        """
    }
    

    process format_association_score {

    publishDir "${output_folder}", overwrite: true, mode: "copy"

    input :
        // Flattening channel so each item can be used independently instead of a bulk
        file REGENIE_PHENO from regenie_output_ch.flatten()
    
    output :
        file "snp_association_*" into snp_asso_ch

    // TODO : uniformize this as a single R script
    """
    # Getting phenotype name
    IN=${REGENIE_PHENO}
    SUFFIX=\${IN##*regenie_output_s2_}
    PHENO=\${SUFFIX%*.regenie}
    
    # Keeping only id and -log10 pvalues + conversion to pvalues for VEGAS compatibility
    Rscript -e "library(magrittr); \\
            readr::read_table('${REGENIE_PHENO}') %>% \\
            dplyr::select(ID, CHISQ) %>% \\
            dplyr::mutate(CHISQ = 10^(-CHISQ)) %>% \\
            readr::write_tsv(file = 'snp_association_\$PHENO', col_names = FALSE)"
    """
    
    }
} else {
    snp_asso_ch = Channel.fromPath("params.snp_association")
}

//////////////////////////////////////////////
///         EXTRACT CONTROLS FOR LD        ///
//////////////////////////////////////////////
if (params.ld_controls == '') {
process extract_controls {

    input:
        file BED from bed
        file bim
        file fam
        path extract_path
        path exclude_path
        path keep_path
        path remove_path

    output:
        file 'plink.bed' into bed_controls
        file 'plink.bim' into bim_controls
        file 'plink.fam' into fam_controls
        
    script:
    """
    plink \\
        --bfile ${BED.baseName} \\
        --filter-controls \\
        --make-bed \\
        $extract \\
        $exclude \\
        $keep \\
        $remove
    """
}
} else {

bed_controls = file("${params.ld_controls}")
bim_controls = file("${bed_controls.getParent()}/${bed_controls.getBaseName()}.bim")
fam_controls = file("${bed_controls.getParent()}/${bed_controls.getBaseName()}.fam")

}

//////////////////////////////////////////////
///                RUN VEGAS               ///
//////////////////////////////////////////////
process vegas {

    tag { "Chr ${CHR}" }

    input:
        file BED from bed_controls
        file bim_controls
        file fam_controls
        file SNPASSOCIATION from snp_asso_ch
        file GLIST from glist_vegas
        val  VEGAS_PARAMS from params.vegas_params
        each CHR from chromosomes

    output:
        file 'scored_genes.vegas.txt' into vegas_out

    """
    vegas2v2 -G -snpandp ${SNPASSOCIATION} -custom `pwd`/${BED.baseName} -glist ${GLIST} -out scored_genes -chr $CHR ${VEGAS_PARAMS}
    sed 's/"//g' scored_genes.out | sed 's/ /\\t/g' >tmp
    R -e 'library(tidyverse); read_tsv("tmp", col_types = "iciddddddcd") %>% filter(!duplicated(Gene)) %>% write_tsv("scored_genes.vegas.txt")'
    """

}

process merge_chromosomes {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file 'scored_genes_chr*' from vegas_out.collect()

    output:
        file 'scored_genes.vegas.txt'

    """
    head -n1 scored_genes_chr1 >scored_genes.vegas.txt
    tail -n +2 -q scored_genes_chr* | sort -n >>scored_genes.vegas.txt
    """

}
