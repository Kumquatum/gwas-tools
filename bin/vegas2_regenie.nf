#!/usr/bin/env nextflow

params.out = '.'

params.gencode = 28
params.genome = '38'
params.vegas_params = ''
params.snp_association = ''
params.ld_controls = ''
params.buffer = 0
params.phenotype = ''
params.covar = ''
params.regenie_params_s1 = ''
params.regenie_params_s2 = ''

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
if (params.snp_association == ''){

  bed = file("${params.bfile}.bed")
  bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
  fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

  if (params.covar == '') {

      process compute_association_score {

          input:
              file BED from bed
              file BIM from bim
              file FAM from fam
              file PHENO from params.phenotype
              val  REGENIE_PARAMS_S1 from params.regenie_params_s1
              val  REGENIE_PARAMS_S2 from params.regenie_params_s2

          output:
              file 'snp_association' into snp_association

          """
          mkdir regenie_run
          regenie \\
            --step 1 \\ 
            --bed ${BED} \\
            --phenoFile ${PHENO} \\
            --bt \\
            --lowmem \\
            --lowmem-prefix regenie_run/tmp_lowmem_regenie_s1 \\
            --out regenie_run/regenie_output_s1 \\
            ${REGENIE_PARAMS_S1}

          regenie \\
            --step 2 \\
            --bed ${BED} \\
            --phenoFile ${PHENO} \\
            --bt \\
            --pred regenie_run/regenie_output_s1_pred.list \\
            --print-pheno \\
            --out regenie_run/regenie_output_s2 \\
            ${REGENIE_PARAMS_S2}
          
          Rscript -e "library(magrittr); \\
            readr::read_table('${REGENIE_OUTPUT_S2}') %>% \\
            dplyr::select(ID, CHISQ) %>% \\
            write_tsv(file = "snp_association", col_names = FALSE)"
          """

      }

  } else {

    process compute_association_score_with_covars {

        input:
            file BED from bed
            file BIM from bim
            file FAM from fam
            file PHENO from params.phenotype
            file COVAR from params.covar
            val  REGENIE_PARAMS_S1 from params.regenie_params_s1
            val  REGENIE_PARAMS_S2 from params.regenie_params_s2

        output:
            file 'snp_association' into snp_association

        """
          mkdir regenie_run
          regenie \\
            --step 1 \\ 
            --bed ${BED} \\
            --phenoFile ${PHENO} \\
            --covarFile ${COVAR} \\
            --bt \\
            --lowmem \\
            --lowmem-prefix regenie_run/tmp_lowmem_regenie_s1 \\
            --out regenie_run/regenie_output_s1 \\
            ${REGENIE_PARAMS_S1}

          regenie \\
            --step 2 \\
            --bed ${BED} \\
            --phenoFile ${PHENO} \\
            --covarFile ${COVAR} \\
            --bt \\
            --pred regenie_run/regenie_output_s1_pred.list \\
            --print-pheno \\
            --out regenie_run/regenie_output_s2 \\
            ${REGENIE_PARAMS_S2}
          
          Rscript -e "library(magrittr); \\
            readr::read_table('${REGENIE_OUTPUT_S2}') %>% \\
            dplyr::select(ID, CHISQ) %>% \\
            write_tsv(file = "snp_association", col_names = FALSE)"        
        """

    }
  }
} else {

  snp_association = file(params.snp_association)

}

//////////////////////////////////////////////
///         EXTRACT CONTROLS FOR LD        ///
//////////////////////////////////////////////
if (params.ld_controls == '') {
  process extract_controls {

      input:
          file BED from bed

      output:
          file 'plink.bed' into bed_controls
          file 'plink.bim' into bim_controls
          file 'plink.fam' into fam_controls

      """
      plink --bfile ${BED.baseName} --filter-controls --make-bed
      """
      // --filter-controls : given case/control data, causes only controls to be included in the current analysis
      // --make-bed        : creates a new PLINK 1 binary fileset, after applying sample/variant filters and other operations 

  }
} else {

  bed_controls = file("${params.ld_controls}.bed")
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
        file SNPASSOCIATION from snp_association
        file GLIST from glist_vegas
        val VEGAS_PARAMS from params.vegas_params
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
