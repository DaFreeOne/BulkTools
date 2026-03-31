FROM rocker/r-ver:4.4.3

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    cmake \
    git \
    pkg-config \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libblas-dev \
    liblapack-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*


# Installing conda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in paths
ENV PATH=$CONDA_DIR/bin:$PATH

# Setting conda env and downloading some R libraries through conda
RUN         "conda create -n 'BulkTools' -c conda-forge r-base=4.4.3 python=3.13.7"
RUN         "conda activate BulkTools"
RUN         "conda install conda-forge::r-shiny conda-forge::r-shinyFiles conda-forge::fs" &&\
            "conda install bioconda::bioconductor-annotationdbi" &&\
            "conda install bioconda::bioconductor-fgsea"

# Installing all required python packages listed in requirements_py.txt
RUN         "pip install -r requirements_py.txt"

# Downloading BiocManager and installing libraries through BiocManager
RUN R -q -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -q -e "BiocManager::install('org.Hs.eg.db', ask=FALSE, update=FALSE)" &&\
    R -q -e "BiocManager::install('tximport', ask=FALSE, update=FALSE)" &&\
    R -q -e "BiocManager::install('GenomeInfoDb', ask=FALSE, update=FALSE)" &&\
    R -q -e "BiocManager::install('AnnotationDbi', ask=FALSE, update=FALSE)" &&\
    R -q -e "BiocManager::install('fgsea', ask=FALSE, update=FALSE)" &&\
    R -q -e "BiocManager::install('GSVA', ask=FALSE, update=FALSE)" &&\
    R -q -e "BiocManager::install('DESeq2', ask=FALSE, update=FALSE)" &&\
    R -q -e "library(org.Hs.eg.db)"

# Installing last R libs
RUN R -q -e "install.packages(c('optparse', 'ggplot2', 'tidyverse', 'dplyr', 'tidyestimate', 'DT', 'bslib', 'jsonlite','data.table', 'Matrix'), repos='https://cloud.r-project.org', Ncpus=max(1, parallel::detectCores()-1))" &&\
    R -q -e "library(optparse); library(tidyestimate); library(jsonlite); library(data.table); library(Matrix)"


# Test si toutes les libs sont bien installées 
RUN R -q -e "library(optparse); library(tidyestimate); library(org.Hs.eg.db); library(jsonlite); library(data.table); library(Matrix); cat('Docker image OK\\n')"

WORKDIR /app
COPY app/ /app/app
COPY REF_DATA/ app/REF_DATA/

EXPOSE 5288

ENTRYPOINT ["Rscript", "/app/demarreur_app.R"]

