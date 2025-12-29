# Use Python 3.10 slim image (small, lightweight)
FROM python:3.10-slim

# Set working directory inside the container
WORKDIR /app

# Copy all files from your project folder into container
COPY . /app

# Install system dependencies needed for bioinformatics packages
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir -r requirements.txt

# Default command to run your script
CMD ["python", "snp_annotation.py"]
