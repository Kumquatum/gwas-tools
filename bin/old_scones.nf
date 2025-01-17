#!/usr/bin/env nextflow

params.out = '.'

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

// SConES parameters
params.network = 'gs'
params.score = 'chi2'
params.criterion = 'stability'
params.encoding = 'additive'
params.eta = null
params.lambda = null

// additional files
snp2gene = (params.network == 'gm' | params.network == 'gi') ? file(params.snp2gene) : file('NO_SNP2GENE')
tab2 = (params.network == 'gi') ? file(params.tab2) : file('NO_TAB2')

process bed2r {

    input:
        file BED from bed
        file BIM from bim
        file FAM from fam

    output:
        file 'gwas.RData' into rgwas_network, rgwas_scones

    script:
    template 'io/bed2r.R'

}

process create_network {

    input:
        file RGWAS from rgwas_network
        val NET from params.network
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        file 'net.RData' into RNET

    """
    #!/usr/bin/env Rscript
    library(martini)
    library(tidyverse)
    library(igraph)

    load("${RGWAS}")
    netType <- "${NET}"

    if (netType == "gs") {
        net <- get_GS_network(gwas)
    } else if (netType %in% c('gm', 'gi')) {
        snp2gene <- read_tsv("${SNP2GENE}")

        if (netType == "gm") {
            net <- get_GM_network(gwas, snpMapping = snp2gene)
        } else if (netType == "gi") {
            tab2 <- read_tsv("${TAB2}", col_types = cols(.default = col_character())) %>%
                rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
                select(gene1, gene2)
            net <- get_GI_network(gwas, snpMapping = snp2gene, ppi = tab2)
        }
    } else {
        stop("network type not recognized.")
    }

    save(net, file = "net.RData")
    """
}

if (params.eta == null | params.lambda == null) {
    process scones {

        publishDir "$params.out", overwrite: true, mode: "copy"

        input:
            file RGWAS from rgwas_scones
            file RNET
            val SCORE from params.score
            val CRITERION from params.criterion

        output:
            file 'cones.tsv'

        script:
        template 'discovery/run_old_scones.R'

        }
} else {
    process parametrized_scones {

        publishDir "$params.out", overwrite: true, mode: "copy"

        input:
            file RGWAS from rgwas_scones
            file RNET
            val SCORE from params.score
            val CRITERION from params.criterion
            val ETA from params.eta
            val LAMBDA from params.lambda

        output:
            file 'cones.tsv'

        script:
        template 'discovery/run_old_scones_params.R'

        }
}
