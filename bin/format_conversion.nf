#!/usr/bin/env nextflow

// Minimal args
params.file_to_convert = ""
params.conversion_type = ""

// Optional args
params.out = "."

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