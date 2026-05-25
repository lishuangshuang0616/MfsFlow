# Installation Guide

This guide provides detailed, step-by-step instructions for installing the MfsFlow toolkit and all its dependencies.

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Quick Installation](#quick-installation)
3. [Installing System Dependencies](#installing-system-dependencies)
4. [Installing Python Dependencies](#installing-python-dependencies)
5. [Installing MfsFlow](#installing-mfsflow)
6. [Verifying Installation](#verifying-installation)
7. [Advanced Installation](#advanced-installation)
8. [Docker Installation](#docker-installation)

## System Requirements

### Minimum Requirements
- **OS**: Linux (Ubuntu 18.04+, CentOS 7+) or macOS (10.14+)
- **CPU**: 4 cores
- **RAM**: 16 GB
- **Storage**: 50 GB free space (for reference indices and temporary files)
- **Python**: 3.8+

### Recommended Requirements
- **CPU**: 20+ cores
- **RAM**: 64+ GB
- **Storage**: 200+ GB SSD
- **Python**: 3.9+

## Quick Installation

For experienced users, here's the quick installation:

```bash
# Clone repository
git clone https://github.com/lishuangshuang0616/MfsFlow.git
cd MfsFlow

# Install Python dependencies
pip install -r requirements.txt

# Install MfsFlow
pip install -e .

# Verify installation
mfsflow --help
```

For detailed instructions, continue reading.

## Installing System Dependencies

**Good news!** MfsFlow comes with bundled bioinformatics tools in the `software/` directory:
- **STAR**: For genome alignment
- **Samtools**: For BAM file manipulation  
- **Pigz**: For parallel gzip compression
- **Seqkit**: For FASTQ manipulation
- **FeatureCounts** (from Subread): For read counting

The pipeline automatically resolves these tools from the `software/` directory on Linux systems. **No manual installation required for most users.**

### Optional: Install System Dependencies Separately

If you prefer to use system-installed versions of these tools, you can install them via Conda or package managers. The pipeline will use system PATH to resolve them.

#### Method 1: Conda (Optional)

```bash
# Create conda environment
conda create -n mfsflow python=3.9
conda activate mfsflow

# Install bioinformatics tools (optional - already bundled)
conda install -c bioconda star samtools subread pigz seqkit

# Verify installations
STAR --version
samtools --version
featureCounts -v
pigz --version
seqkit version
```

#### Method 2: System Package Manager (Optional)

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install samtools pigz
```

**CentOS/RHEL:**
```bash
sudo yum install samtools pigz
```

**macOS:**
```bash
brew install samtools pigz seqkit
```

**Note**: STAR and FeatureCounts may need to be installed manually or via Conda as they might not be available in all system package managers.

## Installing Python Dependencies

### Method 1: pip (Recommended)

```bash
# Create virtual environment (recommended)
python3 -m venv mfsflow_env
source mfsflow_env/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### Method 2: Conda

```bash
# Create conda environment
conda create -n mfsflow python=3.9
conda activate mfsflow

# Install Python dependencies
pip install -r requirements.txt
```

### Python Dependency Details

| Package | Version | Purpose |
|---------|---------|---------|
| PyYAML | >=5.1 | YAML configuration file parsing |
| pysam | >=0.21 | BAM file manipulation |
| numpy | >=1.20 | Numerical computations |
| pandas | >=1.3 | Data manipulation |
| scipy | >=1.7 | Scientific computing |
| matplotlib | >=3.5 | Plotting |
| anndata | >=0.8 | AnnData file format |
| h5py | >=3.1 | HDF5 file support |

**Note**: All packages use open version ranges for maximum compatibility. If you encounter dependency conflicts, you can adjust versions as needed.

### Troubleshooting Python Dependencies

#### Issue: pysam installation fails

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get install zlib1g-dev libcurl4-openssl-dev libssl-dev

# Retry pysam installation
pip install pysam==0.21.0
```

#### Issue: numpy/pandas version conflicts

```bash
# Create fresh environment
python3 -m venv mfsflow_env_clean
source mfsflow_env_clean/bin/activate
pip install -r requirements.txt
```

## Installing MfsFlow

### Method 1: Development Mode (Recommended)

```bash
# Clone repository
git clone https://github.com/lishuangshuang0616/MfsFlow.git
cd MfsFlow

# Install in development mode
pip install -e .

# Verify installation
mfsflow --help
```

### Method 2: Direct Installation

```bash
# Clone repository
git clone https://github.com/lishuangshuang0616/MfsFlow.git
cd MfsFlow

# Install
pip install .

# Verify installation
mfsflow --help
```

### Method 3: From Source (Without pip)

```bash
# Clone repository
git clone https://github.com/lishuangshuang0616/MfsFlow.git
cd MfsFlow

# Add to PYTHONPATH
export PYTHONPATH=/path/to/MfsFlow:$PYTHONPATH

# Use directly
python3 run_analysis_pipeline.py --help
```

## Verifying Installation

### 1. Check MfsFlow CLI

```bash
mfsflow --help
```

Expected output:
```
usage: mfsflow [-h] --fastqs FASTQS --genomeDir GENOMEDIR --sample SAMPLE ...
```

### 2. Check Bundled System Dependencies

MfsFlow automatically uses bundled tools from the `software/` directory. To verify they are accessible:

```bash
# Check STAR (bundled)
ls -lh /path/to/MfsFlow/software/star/bin/STAR

# Check Samtools (bundled)
ls -lh /path/to/MfsFlow/software/samtools/bin/samtools

# Check FeatureCounts (bundled)
ls -lh /path/to/MfsFlow/software/subread/bin/featureCounts

# Check Pigz (bundled)
ls -lh /path/to/MfsFlow/software/pigz/bin/pigz

# Check Seqkit (bundled)
ls -lh /path/to/MfsFlow/software/seqkit/bin/seqkit
```

**Note**: The pipeline automatically resolves these tools. You don't need to add them to your PATH.

To verify the pipeline can find them, run a test analysis (see Step 4).

### 3. Check Python Dependencies

```bash
# Run test suite
python3 -m unittest discover -s tests
```

All tests should pass.

### 4. Run a Test Analysis

```bash
# Download test data (if available)
# ... instructions for test data ...

# Run pipeline on test data
mfsflow \
  --fastqs /path/to/test/fastqs \
  --genomeDir /path/to/test/reference \
  --sample test_sample \
  --plate 1 \
  --outdir /path/to/test/output \
  --threads 4
```

## Advanced Installation

### Installing to Custom Location

```bash
# Install to custom prefix
pip install --prefix=/custom/install/path .

# Add to PATH
export PATH=/custom/install/path/bin:$PATH
export PYTHONPATH=/custom/install/path/lib/python3.9/site-packages:$PYTHONPATH
```

### Installing Specific Version

```bash
# Clone specific version
git clone https://github.com/lishuangshuang0616/MfsFlow.git
cd MfsFlow
git checkout v0.1.0

# Install
pip install -e .
```

### Offline Installation

```bash
# Download dependencies on internet-connected machine
pip download -r requirements.txt -d /path/to/wheelhouse

# Transfer wheelhouse to target machine

# Install offline
pip install --no-index --find-links=/path/to/wheelhouse -r requirements.txt
pip install --no-index --find-links=/path/to/wheelhouse .
```

## Docker Installation

### Using Pre-built Image (If Available)

```bash
# Pull Docker image
docker pull mfsflow/mfsflow:latest

# Run container
docker run -v /path/to/data:/data mfsflow/mfsflow:latest \
  mfsflow --fastqs /data/fastqs --genomeDir /data/reference ...
```

### Building Docker Image

Create `Dockerfile`:

```dockerfile
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    gcc \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools
RUN apt-get update && apt-get install -y \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz \
    && tar -xzf 2.7.10b.tar.gz \
    && cp STAR-2.7.10b/source/STAR /usr/local/bin/ \
    && rm -rf 2.7.10b.tar.gz STAR-2.7.10b

# Install FeatureCounts
RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-source.tar.gz \
    && tar -xzf subread-2.0.3-source.tar.gz \
    && cd subread-2.0.3-source/src \
    && make -f Makefile.Linux \
    && cp ../bin/featureCounts /usr/local/bin/ \
    && rm -rf /subread-2.0.3-source*

# Install Pigz and Seqkit
RUN apt-get update && apt-get install -y pigz && rm -rf /var/lib/apt/lists/*
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.3.0/seqkit_linux_amd64.tar.gz \
    && tar -xzf seqkit_linux_amd64.tar.gz \
    && mv seqkit /usr/local/bin/ \
    && rm seqkit_linux_amd64.tar.gz

# Copy MfsFlow
COPY . /mfsflow
WORKDIR /mfsflow

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir .

# Set entrypoint
ENTRYPOINT ["mfsflow"]
```

Build and run:

```bash
# Build image
docker build -t mfsflow:latest .

# Run container
docker run -v /path/to/data:/data mfsflow:latest \
  --fastqs /data/fastqs \
  --genomeDir /data/reference \
  --sample SAMPLE \
  --plate 1 \
  --outdir /data/output \
  --threads 20
```

## Next Steps

After installation:

1. Proceed to [User Guide](user_guide.md) for usage instructions
2. Prepare reference genome (see User Guide)
3. Run test analysis to verify installation

## Getting Help

If you encounter installation issues:

1. Check [Troubleshooting Guide](troubleshooting.md)
2. Open an issue on GitHub
3. Contact the development team
