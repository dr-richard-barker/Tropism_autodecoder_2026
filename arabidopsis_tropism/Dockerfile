FROM condaforge/mambaforge:latest

LABEL maintainer="Phylo Biomni"
LABEL description="Arabidopsis thaliana spaceflight tropism recognition pipeline"
LABEL version="1.0"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    graphviz \
    libgraphviz-dev \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment
COPY environment.yml /opt/arabidopsis-tropism/environment.yml
WORKDIR /opt/arabidopsis-tropism
RUN mamba env create -f environment.yml && mamba clean -a

# Install ggPlantmap from GitHub
RUN conda run -n arabidopsis-tropism Rscript -e \
    'remotes::install_github("leonardojo/ggPlantmap", upgrade="never")'

# Copy code
COPY Code/ /opt/arabidopsis-tropism/Code/
COPY Data/ /opt/arabidopsis-tropism/Data/

# Activate environment
ENV PATH /opt/conda/envs/arabidopsis-tropism/bin:$PATH

# Entry point
ENTRYPOINT ["python", "/opt/arabidopsis-tropism/Code/train_autodecoder.py"]
CMD ["--help"]
