# =============================================================================
# Dockerfile — OmicsFlow Multi-Omics Pipeline
# =============================================================================
#
# BUILD:
#   docker build -t omicsflow:2.0.1 .
#
# RUN (standalone):
#   docker run --rm -v $(pwd)/data:/omicsflow/data -v $(pwd)/results:/omicsflow/results \
#     omicsflow:2.0.1 Rscript modules/rna/preprocess_rna.R --input data/rna_expression_raw.rds
#
# RUN (via Nextflow):
#   nextflow run main.nf -profile docker --rna data/rna_expression_raw.rds ...
#
# =============================================================================

# --- Stage 1: Base R + System Libraries ---
FROM rocker/r-ver:4.4.2

LABEL maintainer="Abdullah Ibrahim Ali"
LABEL version="2.0.1"
LABEL description="OmicsFlow: Multi-Omics Integration & Discovery Pipeline"

# Prevent interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# --- System Dependencies ---
# Required for: R package compilation, HDF5 (MOFA2), Python (mofapy2),
# graphics (Cairo/png), SSL/curl, XML parsing, git
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    # Compilation toolchain
    build-essential \
    cmake \
    pkg-config \
    gfortran \
    # Networking & crypto
    libcurl4-openssl-dev \
    libssl-dev \
    # XML & text
    libxml2-dev \
    libpcre2-dev \
    # HDF5 (required for MOFA2 model I/O)
    libhdf5-dev \
    # Graphics rendering
    libcairo2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    # Compression
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    # Git (for remotes::install_github if needed)
    libgit2-dev \
    git \
    # Python 3 (required for MOFA2 mofapy2 backend)
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    # Misc
    procps \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# --- Python: Install mofapy2 (MOFA2 backend) ---
RUN pip3 install --no-cache-dir --break-system-packages \
    mofapy2==0.7.1 \
    numpy \
    scipy \
    pandas \
    scikit-learn \
    h5py

# =============================================================================
# Stage 2: R Package Installation
# =============================================================================

# --- BiocManager (pinned Bioconductor version) ---
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(version = '3.20', ask = FALSE, update = FALSE)"

# --- CRAN Packages (Layer 1: Core utilities) ---
RUN R -e "install.packages(c( \
    'optparse', \
    'jsonlite', \
    'data.table', \
    'matrixStats', \
    'stringr', \
    'dplyr', \
    'tidyr', \
    'tibble', \
    'scales', \
    'RColorBrewer', \
    'pheatmap' \
), repos='https://cloud.r-project.org', Ncpus=4)"

# --- CRAN Packages (Layer 2: Visualization) ---
RUN R -e "install.packages(c( \
    'ggplot2', \
    'ggrepel', \
    'cowplot', \
    'gridExtra' \
), repos='https://cloud.r-project.org', Ncpus=4)"

# --- CRAN Packages (Layer 3: ML & Survival) ---
RUN R -e "install.packages(c( \
    'survival', \
    'caret', \
    'glmnet', \
    'xgboost', \
    'randomForestSRC' \
), repos='https://cloud.r-project.org', Ncpus=4)"

# --- Bioconductor Packages (Layer 4: Core genomics) ---
RUN R -e "BiocManager::install(c( \
    'SummarizedExperiment', \
    'edgeR', \
    'limma', \
    'sva', \
    'impute' \
), ask = FALSE, update = FALSE, Ncpus=4)"

# --- Bioconductor Packages (Layer 5: Methylation annotation) ---
RUN R -e "BiocManager::install(c( \
    'IlluminaHumanMethylation450kanno.ilmn12.hg19', \
    'minfi' \
), ask = FALSE, update = FALSE, Ncpus=4)"

# --- Bioconductor Packages (Layer 6: Genomic ranges) ---
RUN R -e "BiocManager::install(c( \
    'GenomicRanges', \
    'IRanges', \
    'GenomeInfoDb', \
    'biomaRt' \
), ask = FALSE, update = FALSE, Ncpus=4)"

# --- Bioconductor Packages (Layer 7: Mutations) ---
RUN R -e "BiocManager::install(c( \
    'maftools' \
), ask = FALSE, update = FALSE, Ncpus=4)"

# --- Bioconductor Packages (Layer 8: Integration — MOFA2) ---
RUN R -e "BiocManager::install(c( \
    'MOFA2', \
    'basilisk', \
    'reticulate', \
    'rhdf5' \
), ask = FALSE, update = FALSE, Ncpus=4)"

# --- Bioconductor Packages (Layer 9: Enrichment) ---
RUN R -e "BiocManager::install(c( \
    'clusterProfiler', \
    'org.Hs.eg.db', \
    'enrichplot', \
    'DOSE', \
    'AnnotationDbi' \
), ask = FALSE, update = FALSE, Ncpus=4)"

# =============================================================================
# Stage 3: Project Setup
# =============================================================================

# Create project structure
RUN mkdir -p /omicsflow

WORKDIR /omicsflow

# Copy pipeline code (R modules, Nextflow, configs)
COPY modules/         /omicsflow/modules/
COPY modules_nf/      /omicsflow/modules_nf/
COPY configs/         /omicsflow/configs/
COPY tests/           /omicsflow/tests/
COPY assets/          /omicsflow/assets/
COPY envs/            /omicsflow/envs/
COPY main.nf          /omicsflow/main.nf
COPY nextflow.config  /omicsflow/nextflow.config

# Create data and results mount points
RUN mkdir -p /omicsflow/data /omicsflow/results

# =============================================================================
# Stage 4: Smoke Test & Metadata
# =============================================================================

# Copy and run smoke test to validate the image
COPY docker/smoke_test.R /omicsflow/docker/smoke_test.R
RUN Rscript /omicsflow/docker/smoke_test.R

# Runtime metadata
ENV OMICSFLOW_version="2.0.1"
ENV BIOCONDUCTOR_VERSION="3.20"
ENV R_VERSION="4.4.2"

# Default entrypoint
ENTRYPOINT ["Rscript"]
CMD ["--help"]
