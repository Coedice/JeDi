# Use mambaorg/micromamba as base for efficient conda package management
FROM mambaorg/micromamba:1.5.10

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Switch to micromamba user
USER $MAMBA_USER

# Install bioinformatics tools via conda
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.13 \
    samtools=1.20 \
    bcftools=1.20 \
    vcftools=0.1.16 \
    stacks \
    bedtools=2.31 \
    mawk=1.3.4 \
    graphviz \
    piawka \
    && micromamba clean --all --yes

# Install uv
USER root
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:$PATH"

# Set working directory
WORKDIR /workspace

# Copy project files
COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace/

# Install Python dependencies with uv
USER root
RUN uv sync

# Default command - stay as root for runtime
CMD ["/bin/bash"]
