# LAMBclassifier in a reproducible container.
#
# Build (from the package root):
#   docker build -t lambclassifier .
# Run on a file or folder in the current directory (results are written next to it):
#   docker run --rm -v "$PWD":/data lambclassifier /data/my_expression.csv
# Try the shipped example (results go to /data/example_out):
#   docker run --rm -v "$PWD":/data --entrypoint Rscript lambclassifier -e \
#     'LAMBclassifier::lamb_run(system.file("extdata","example_bulk.csv",package="LAMBclassifier"), out_dir = "/data/example_out")'
#
# Pinned: R 4.4.2 (rocker/r-ver installs CRAN packages from a dated Posit
# Package Manager snapshot), BiocManager 1.30.25, Bioconductor 3.20.

FROM rocker/r-ver:4.4.2

ARG BIOC_VERSION=3.20
ARG BIOCMANAGER_VERSION=1.30.25
# Optional readers named in Suggests (for example "SeuratObject SingleCellExperiment").
ARG EXTRA_PACKAGES=""

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev libhdf5-dev \
        libpng-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages('remotes')" \
 && Rscript -e "remotes::install_version('BiocManager', version = '${BIOCMANAGER_VERSION}', upgrade = 'never')" \
 && Rscript -e "BiocManager::install(version = '${BIOC_VERSION}', ask = FALSE, update = FALSE)" \
 && Rscript -e "BiocManager::install(c('singscore', 'genefu'), version = '${BIOC_VERSION}', ask = FALSE, update = FALSE); \
                stopifnot(requireNamespace('singscore'), requireNamespace('genefu'))"

COPY . /tmp/LAMBclassifier
RUN Rscript -e "options(repos = BiocManager::repositories()); \
                remotes::install_local('/tmp/LAMBclassifier', dependencies = NA, upgrade = 'never'); \
                extra <- strsplit('${EXTRA_PACKAGES}', ' ')[[1]]; extra <- extra[nzchar(extra)]; \
                if (length(extra)) BiocManager::install(extra, ask = FALSE, update = FALSE); \
                stopifnot(requireNamespace('LAMBclassifier')); \
                print(LAMBclassifier::lamb_specification()[c('version', 'signature_sha256')])" \
 && rm -rf /tmp/LAMBclassifier

WORKDIR /data
ENTRYPOINT ["Rscript", "-e", "a <- commandArgs(TRUE); if (!length(a)) stop('Give a file or folder, e.g. docker run --rm -v \"$PWD\":/data lambclassifier /data/my_data.csv', call. = FALSE); LAMBclassifier::lamb_run(a[1])"]
