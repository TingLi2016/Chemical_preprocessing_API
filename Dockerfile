# ---- Base Image with Conda ----
FROM continuumio/miniconda3:latest AS base

# Install curl for healthchecks
RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

ENV CONDA_ENV_NAME=chem-preprocessing \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

# Copy environment and project metadata
COPY environment.yml pyproject.toml ./ 
COPY . .

# Create environment
RUN conda update -n base -c defaults conda -y && \
    conda env create -n ${CONDA_ENV_NAME} -f environment.yml && \
    conda clean -afy

# Use conda env for all RUN commands from here
SHELL ["conda", "run", "-n", "chem-preprocessing", "/bin/bash", "-c"]

# Install editable project and dev dependencies
RUN pip install -e .[all]

# ---- Backend Image ----
FROM base AS backend

EXPOSE 8008
CMD ["conda", "run", "-n", "chem-preprocessing", "uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8008", "--log-level", "info"]

# ---- Frontend Image with Python 3.10 ----
FROM python:3.10-slim AS frontend

WORKDIR /app

# Install system dependencies for Streamlit and build tools
RUN apt-get update && apt-get install -y \
    build-essential \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    curl && \
    rm -rf /var/lib/apt/lists/*

# Copy only what’s needed
COPY pyproject.toml environment.yml ./ 
COPY . .

# Install Python dependencies
RUN pip install --upgrade pip && \
    pip install -e .[frontend]

EXPOSE 8502
CMD ["streamlit", "run", "frontend.py", "--server.port=8502", "--server.address=0.0.0.0"]
