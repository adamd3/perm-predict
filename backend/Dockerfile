# Stage 1: mamba_builder - Install alvadesccliwrapper and alvaDescCLI
FROM mambaorg/micromamba:1.5.1 AS mamba_builder

# Set up micromamba environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /tmp/alvadesc_build

# Create mamba environment and install dependencies
COPY alvaDescCLIWrapper ./alvaDescCLIWrapper
RUN micromamba create -n alvadesc python=3.8 openjdk=11 -c conda-forge -c bioconda -c defaults && \
    micromamba run -n alvadesc bash -c " \
    cp -R alvaDescCLIWrapper/alvadesccliwrapper /opt/conda/envs/alvadesc/lib/python3.8/site-packages/ \
    " && \
    micromamba clean --all --yes

# Stage 2: backend_runner - Build the final backend image
FROM python:3.8-slim-bullseye AS backend_runner

WORKDIR /app

# Install system dependencies for RDKit and curl
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy alvaDescCLI executable from the build context (project root)
COPY alvaDescCLI_linux /app/backend/app/ml_models/alvaDescCLI
COPY ./backend/test_feature_generation.py /app/test_feature_generation.py

# Copy alvadesccliwrapper package from mamba_builder stage
# This assumes alvadesccliwrapper is installed in site-packages
# We need to find the exact path in the mamba environment
# For now, let's assume a common path, but this might need adjustment
COPY --from=mamba_builder /opt/conda/envs/alvadesc/lib/python3.8/site-packages/alvadesccliwrapper /usr/local/lib/python3.8/site-packages/alvadesccliwrapper
COPY --from=mamba_builder /opt/conda/envs/alvadesc/lib/python3.8/site-packages/alvadesccliwrapper-*.egg-info /usr/local/lib/python3.8/site-packages/

COPY ./backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY ./backend/app ./app
COPY alvadesc_descriptor_names.json /app/alvadesc_descriptor_names.json
COPY full_feature_names.json /app/full_feature_names.json
COPY desc_smiles.tsv /app/desc_smiles.tsv

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
