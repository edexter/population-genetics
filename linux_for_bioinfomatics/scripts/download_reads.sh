#!/bin/bash

# Download and prepare E. coli sequencing reads for Exercise 1
# This script downloads one sample and downsamples it for faster processing

# ============================================================================
# IMPORTANT: Set your working directory
# ============================================================================
# Change this path to match YOUR project folder location
WORKING_DIR="/home/eric/projects/population-genetics"

# ============================================================================
# Configuration
# ============================================================================
SRA_ACCESSION="SRR11816069"  # E. coli STEC sample from PRJNA633966
TARGET_READS=100000          # Number of read pairs (affects file size and coverage)
SEED=100                     # Random seed for reproducibility

# ============================================================================
# Setup directories
# ============================================================================
echo "Setting up directory structure..."
mkdir -p ${WORKING_DIR}/exercise1/data/reads/downsampled
cd ${WORKING_DIR}/exercise1/data/reads

# ============================================================================
# Check required software
# ============================================================================
echo "Checking for required software..."

command -v wget >/dev/null 2>&1 || {
    echo "ERROR: wget not found. Please install wget."
    exit 1
}

command -v seqtk >/dev/null 2>&1 || {
    echo "ERROR: seqtk not found. Install with:"
    echo "  conda install -c bioconda seqtk"
    exit 1
}

echo "All required software found."
echo ""

# ============================================================================
# Get download URLs from ENA
# ============================================================================
echo "=========================================="
echo "Downloading E. coli sequencing reads"
echo "Sample: ${SRA_ACCESSION}"
echo "=========================================="
echo ""

echo "Getting download URLs from ENA..."
URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRA_ACCESSION}&result=read_run&fields=fastq_ftp" | tail -1)

if [ -z "$URLS" ] || [ "$URLS" == "fastq_ftp" ]; then
    echo "ERROR: Could not get URLs for ${SRA_ACCESSION}"
    exit 1
fi

URL1="ftp://$(echo $URLS | cut -d';' -f1)"
URL2="ftp://$(echo $URLS | cut -d';' -f2)"

# ============================================================================
# Download reads
# ============================================================================
echo "Downloading reads from ENA (European Nucleotide Archive)..."
echo "This may take a few minutes..."
echo ""

wget -O ${SRA_ACCESSION}_1.fastq.gz "$URL1"
wget -O ${SRA_ACCESSION}_2.fastq.gz "$URL2"

# Check if download succeeded
if [ ! -f "${SRA_ACCESSION}_1.fastq.gz" ] || [ ! -s "${SRA_ACCESSION}_1.fastq.gz" ]; then
    echo "ERROR: Download failed"
    exit 1
fi

echo ""
echo "Download complete!"

# ============================================================================
# Downsample reads
# ============================================================================
echo ""
echo "Downsampling to ${TARGET_READS} read pairs..."

seqtk sample -s${SEED} ${SRA_ACCESSION}_1.fastq.gz ${TARGET_READS} | gzip > downsampled/sample_R1.fastq.gz
seqtk sample -s${SEED} ${SRA_ACCESSION}_2.fastq.gz ${TARGET_READS} | gzip > downsampled/sample_R2.fastq.gz

# Remove full files
rm ${SRA_ACCESSION}_1.fastq.gz ${SRA_ACCESSION}_2.fastq.gz

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=========================================="
echo "Download and downsampling complete!"
echo "=========================================="

READS=$(zcat downsampled/sample_R1.fastq.gz | wc -l | awk '{print $1/4}')
SIZE_R1=$(du -h downsampled/sample_R1.fastq.gz | cut -f1)
SIZE_R2=$(du -h downsampled/sample_R2.fastq.gz | cut -f1)

echo "Sample: ${SRA_ACCESSION}"
echo "Read pairs: ${READS}"
echo "File sizes: ${SIZE_R1} (R1), ${SIZE_R2} (R2)"
echo ""
echo "Files saved to:"
echo "  ${WORKING_DIR}/exercise1/data/reads/downsampled/sample_R1.fastq.gz"
echo "  ${WORKING_DIR}/exercise1/data/reads/downsampled/sample_R2.fastq.gz"
echo ""
echo "You can now proceed with read alignment in the next step."
echo "=========================================="
