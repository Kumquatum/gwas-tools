########################
# Setting up basic env #
########################
FROM r-base

RUN apt-get update \
    && apt-get install -y wget unzip tar python2 python3 python3-pip  sed bedtools libcurl4-gnutls-dev libxml2-dev libssl-dev gawk \
    && rm -rf /var/lib/apt/lists/*
RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python2 get-pip.py \
    && rm get-pip.py
RUN mkdir /gwas-tools
ENV PATH="/gwas-tools:${PATH}"
WORKDIR /gwas-tools

#######################################
# Installing data preprocessing tools #
#######################################

# TODO : rework the R package install to get comon dependencies on one side, and each network tool and its dependencies with its own install line + comment
# R packages useful for pipelines
RUN R -e "install.packages(c('mvtnorm', 'corpcor', 'tidyverse', 'magrittr', 'BiocManager', 'cowplot'), repos = 'http://cran.us.r-project.org')"
# VEGAS2 (by impersonating a browser bc not the right file downloaded if done without this)
RUN wget --user-agent="Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1)" https://vegas2.qimrberghofer.edu.au/vegas2v2 \ 
    && chmod a+x vegas2v2
# PLINK
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip \
    && unzip plink_linux_x86_64_20200616.zip
# Regenie 
RUN wget https://github.com/rgcgithub/regenie/releases/download/v3.1.2/regenie_v3.1.2.gz_x86_64_Linux_mkl.zip \
    && unzip regenie_v3.1.2.gz_x86_64_Linux_mkl.zip \
    && rm regenie_v3.1.2.gz_x86_64_Linux_mkl.zip \
    && mv regenie_v3.1.2.gz_x86_64_Linux_mkl regenie

########################################
# Installing network-guided gwas tools #
########################################

# Network tool - SConES (martini) and its dependencies
RUN R -e "BiocManager::install(c('martini','BioNet','twilight'))"
# Network tool - SigMod
RUN wget https://github.com/YuanlongLiu/SigMod/blob/20c561876d87a0faca632a6b93882fcffd719b17/SigMod_v2.zip \
    && unzip SigMod_v2.zip
# Network tool - HotNet2
RUN wget https://github.com/raphael-group/hotnet2/archive/refs/tags/v1.2.1.tar.gz \
    && tar xzf v1.2.1.tar.gz \
    && rm v1.2.1.tar.gz \
    && mv hotnet2-1.2.1 hotnet2 \
    && pip2 install -r hotnet2/requirements.txt
# Network tool - dmGWAS
RUN wget --no-check-certificate https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0.tar.gz \
    && R -e 'install.packages("dmGWAS_3.0.tar.gz", repos = NULL, type="source")'
# Other R packages
RUN R -e "install.packages(c('ranger','SKAT','biglasso','bigmemory','igraph','LEANR','CASMAP', 'doMC'), repos = 'http://cran.us.r-project.org')"


#######################################################
# Installing network-guided epistasis detection tools #
#######################################################

# BEDOPS
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 \
    && tar jxvf bedops_linux_x86_64-v2.4.35.tar.bz2 \
    && cp bin/* .


# Finalizing
WORKDIR /home/docker
