# Use the official Python image as the base image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Install system dependencies required for certain Python packages (e.g., RDKit, TensorFlow, etc.)
RUN apt-get update && apt-get install -y \
    libarchive-dev \
    build-essential \
    libjpeg-dev \
    libpng-dev \
    libxrender-dev \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip to ensure we're using the latest version
RUN pip install --upgrade pip

# Install DeepChem and its dependencies via pip
RUN pip install deepchem==2.8.0 \
    tensorflow==2.12.0 \
    scikit-learn==1.2.2 \
    rdkit==2023.3.3 \
    numpy==1.23 \
    pandas==2.2.3 \
    matplotlib==3.9.2 \
    optuna==4.1.0 \
    networkx==3.2.1

# Install DGL version 1.1.2+cu118 with pip
RUN pip install dgl==1.1.2+cu118 -f https://data.dgl.ai/wheels/cu118/repo.html

# Install PyTorch with CUDA 11.8 support via pip
RUN pip install torch==2.5.0 torchvision==0.20.0 torchaudio==2.5.0 --index-url https://download.pytorch.org/whl/cu118

RUN pip install pyg-lib -f https://data.pyg.org/whl/torch-2.5.0+cu118.html

# Install additional dependencies via pip (PyG, RDKit, etc.)
RUN pip install torch-scatter torch-sparse torch-cluster torch-spline-conv \
    torch_geometric rdkit ipykernel dgllife

# Copy the current directory contents into the container
COPY . /app

# Set the default command to run the application
CMD ["python"]




