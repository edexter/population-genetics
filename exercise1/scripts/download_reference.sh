#!/bin/bash

# Download and prepare E. coli K-12 MG1655 reference genome
# This script downloads the reference genome from NCBI and prepares it for analysis

# ============================================================================
# IMPORTANT: Set your working directory
# ============================================================================
# Change this path to match YOUR project folder location
WORKING_DIR="/home/eric/projects/population-genetics"

# ============================================================================
# Setup directories
# ============================================================================
echo "Setting up directory structure..."
mkdir -p ${WORKING_DIR}/exercise1/data/genomes
cd ${WORKING_DIR}/exercise1/data/genomes

# ============================================================================
# Download reference genome
# ============================================================================
echo "Downloading E. coli K-12 MG1655 reference genome from NCBI..."
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

# ============================================================================
# Decompress
# ============================================================================
echo "Decompressing genome file..."
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

# ============================================================================
# Rename for simplicity
# ============================================================================
echo "Renaming to ecoli_K12_MG1655.fasta..."
mv GCF_000005845.2_ASM584v2_genomic.fna ecoli_K12_MG1655.fasta

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "============================================"
echo "Reference genome download complete!"
echo "============================================"
echo "Genome location: ${WORKING_DIR}/exercise1/data/genomes/ecoli_K12_MG1655.fasta"
echo "Genome size:"
wc -l ecoli_K12_MG1655.fasta
echo ""
echo "You can now use this reference for alignment in the next steps."
echo "============================================"
