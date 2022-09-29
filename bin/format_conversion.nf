#!/usr/bin/env nextflow

// Minimal args
params.file_to_convert = ""
params.conversion_type = ""

// Optional args
params.out = "."
params.additional_file = ""

// 
conv_type = params.conversion_type

if (params.conversion_type == "bed2gen") {

    bed = file("${params.file_to_convert}")
    bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
    fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

    process bed2gen {

        input:
            file BED from bed
            file BIM from bim // Not explicitely used be it is
            file FAM from fam // Idem
        output:
            file 'out.gen' into gen_ch
            file 'out.sample' into sample_ch
        
        script:
            template 'io/bed2gen.sh'
    }

    process rename {

        publishDir "$params.out", overwrite: true, mode: "copy"

        input:
            file BED from bed
            file GEN from gen_ch
            file SAMPLE from sample_ch
        
        output:
            file '*.gen'
            file '*.sample'

        script:
            """
            mv out.gen ${BED.baseName}.gen
            mv out.sample ${BED.baseName}.sample
            """
    }
} else if (params.conversion_type == "fam2phenotype") {

    fam = file("${params.file_to_convert}")

    process fam2phenotype {
        publishDir "$params.out", overwrite: true, mode: "copy"

        input:
            file FAM from fam
            val PHENO from fam.baseName

        output:
            file '*.tsv'
        
        script:
            template 'io/fam2phenotype.sh'
    }
} else if (params.conversion_type == "vegas2hgnc" | params.conversion_type == "vegas2ensembl") {

    vegas = file("${params.file_to_convert}")
    ref = file("${params.additional_file}")
    // Keeping only the id reference format. Either 'hgnc' or 'ensembl'
    converter = params.conversion_type.replaceFirst(/vegas2/, "")

    process vegas2other_gene_ref {
        echo true
        publishDir "$params.out", overwrite: true, mode: "copy"

        input:
            file VEGAS from vegas
            file REF from ref
            val CONV_TO from converter 
        
        output:
            file '*.txt'
        
        script:
        // Retrieving the id reference format in the file before conversion by
        // oppposition to the one provided with CONV_TO
        CONV_FROM = CONV_TO == "hgnc" ? "ensembl" : "hgnc"
            """
            #!/usr/bin/env Rscript
            library(magrittr)

            vegas <- readr::read_table("${VEGAS}")
            ref <- readr::read_table("${REF}")

            merge <- vegas %>%
                dplyr::left_join(ref, by = c("Gene" = "${CONV_FROM}_gene_id")) %>%
                dplyr::mutate(Gene = ${CONV_TO}_gene_id) %>%
                dplyr::select(-snp, -${CONV_TO}_gene_id)

            if (any(is.na(merge\$Gene))) 
                # Using cat instead of warning as stderr can't be outputed in stdout via `echo true` in nextflow
                # warning("Some ${CONV_FROM} IDs didn't mached ${CONV_TO} IDs. Removing them")
                cat("\\nSome ${CONV_FROM} IDs didn't mached ${CONV_TO} IDs. Removing them")

            merge %>%
                dplyr::distinct() %>%
                dplyr::filter(!is.na(Gene)) %>%
                readr::write_tsv("${VEGAS.baseName}_conv_${CONV_TO}.txt")
            """
            // TODO: Find out how to display the warning since `debug true` isn't working
    }
} else {
    process no_conversion {
        input:
            var CONV_TYPE from conv_type
        
        """
        echo 'Unknown conversion_type: ${CONV_TYPE}'
        exit 1
        """
    }
}