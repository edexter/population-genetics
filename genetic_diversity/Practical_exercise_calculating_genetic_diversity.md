# Practical Exercise: Calculating Genetic Diversity from VCF Data

## Learning Objectives

By the end of this practical, you will be able to:
1. Set up a proper bioinformatics project directory structure
2. Inspect and filter VCF files containing genomic variant data
3. Calculate two measures of genetic diversity: S (segregating sites) and π (nucleotide diversity)
4. Interpret these statistics in the context of population genetics theory
5. Visualize diversity patterns across a genome

---

## Background

In the lecture, we learned about several statistics for measuring genetic diversity. In this practical, we'll focus on calculating two fundamental measures using bioinformatics tools:

- **S**: Number of segregating (polymorphic) sites - the simplest measure of variation
- **π (pi)**: Nucleotide diversity - average pairwise differences between sequences

We'll be working with a VCF (Variant Call Format) file containing SNP data from multiple samples of the bacterium *Pasteuria ramosa* collected from a pond population.

**Note about our organism:** *Pasteuria ramosa* is a haploid bacterium, meaning each cell has only one copy of each gene. When we calculate diversity statistics, we're measuring variation across the population, not within individuals.

---

## Section 1: Workspace Preparation

### 1.1 Check Required Software

Before we begin, we need to ensure you have the necessary bioinformatics tool installed. We'll be using:
- **vcftools**: For VCF filtering and calculating diversity statistics (required)
- **bcftools**: Optional tool with sometimes more convenient features for VCF manipulation

#### Check if vcftools is installed:

```bash
# Check vcftools
vcftools --version
```

**Expected output:**
- vcftools: version 0.1.16 or higher

#### If vcftools is not installed:

**On Ubuntu/Debian Linux:**
```bash
sudo apt-get update
sudo apt-get install vcftools
```

**On macOS with Homebrew:**
```bash
brew install vcftools
```

**On a computing cluster:**
Check with your system administrator or use module commands:
```bash
module avail vcftools
module load vcftools
```

#### Optional: Install bcftools

Throughout this tutorial, we'll show alternative commands using bcftools where applicable. While not required, bcftools can sometimes offer more concise syntax for certain operations.

**Installation (same as vcftools):**

**On Ubuntu/Debian Linux:**
```bash
sudo apt-get install bcftools
```

**On macOS with Homebrew:**
```bash
brew install bcftools
```

**On a computing cluster:**
```bash
module load bcftools
```

---

### 1.2 Create Project Directory Structure

A well-organized project directory is essential for reproducible bioinformatics work. We'll create a standard structure:

```
genetic_diversity/
├── data/           # Raw data files (VCF)
├── scripts/        # Analysis scripts
├── results/        # Output files and statistics
└── scratch/        # Temporary files
```

#### Create the directory structure:

```bash
# Navigate to your working directory (modify as needed)
cd ~

# Create main project directory
mkdir -p genetic_diversity

# Navigate into it
cd genetic_diversity

# Create subdirectories
mkdir -p data scripts results scratch

# Verify the structure
ls -l
```

**Expected output format:**
```
drwxr-xr-x  2 user group 4096 Oct 30 10:00 data
drwxr-xr-x  2 user group 4096 Oct 30 10:00 results
drwxr-xr-x  2 user group 4096 Oct 30 10:00 scratch
drwxr-xr-x  2 user group 4096 Oct 30 10:00 scripts
```

---

## Section 2: Obtain and Organize the Data

### 2.1 Download the VCF File

The VCF file contains variant information for multiple *Pasteuria ramosa* samples from a single pond population. This file has been pre-processed and filtered for quality.

**Download the data:**

```bash
# Make sure you're in the project directory
cd ~/genetic_diversity

# Download the VCF file from GitHub using wget
wget -O data/pasteuria_population.vcf https://raw.githubusercontent.com/edexter/population-genetics/main/genetic_diversity/data/pasteuria_population.vcf

# Alternative: Download using curl
# curl -o data/pasteuria_population.vcf https://raw.githubusercontent.com/edexter/population-genetics/main/genetic_diversity/data/pasteuria_population.vcf
```

**Note:** We're downloading the raw file from GitHub, which is why we use the `raw.githubusercontent.com` URL.

---

### 2.2 Verify the File

Let's make sure the file downloaded correctly and is in the right location:

```bash
# Check the file exists and see its size
ls -lh data/pasteuria_population.vcf

# View the first few lines
head -n 20 data/pasteuria_population.vcf
```

**What you should see:**
- VCF files start with header lines beginning with `##`
- The last header line starts with `#CHROM` and lists all sample names
- Data lines follow, one per variant site

---

## Section 3: Explore and Validate the VCF

Before calculating diversity statistics, we need to understand what's in our VCF file. This step is crucial for quality control.

### 3.1 Basic VCF Statistics

Let's get an overview of the data using vcftools:

```bash
# Get basic statistics about the VCF file
vcftools --vcf data/pasteuria_population.vcf

# This will output:
# - Number of individuals (samples)
# - Number of sites (variants) after filtering
```

**Questions to answer:**
- How many variant sites are in the file?
- How many samples are included?

---

### 3.2 Detailed Statistics

We can look at the VCF file header to see sample names and understand the data format:

```bash
# View VCF header including sample names
grep "^#CHROM" data/pasteuria_population.vcf

# Count how many samples (columns after the first 9 fixed columns)
grep "^#CHROM" data/pasteuria_population.vcf | awk '{print "Number of samples:", NF-9}'

# Look at the first few data lines
grep -v "^##" data/pasteuria_population.vcf | head -n 5
```

---

### 3.3 Filter to Remove Indels and Keep Only Biallelic SNPs

Most population genetic statistics are designed for SNPs (single nucleotide polymorphisms) only. We need to:
1. Remove indels (insertions/deletions)
2. Keep only biallelic sites (exactly 2 alleles: reference + 1 alternate)

We'll use a systematic naming convention to track our filtering:
- `pasteuria_population.vcf` → original file
- `pasteuria_population_filtered.vcf` → filtered file (SNPs only, biallelic)

**Method 1: Using vcftools (recommended for this course)**

```bash
# Filter to keep only biallelic SNPs using vcftools
vcftools --vcf data/pasteuria_population.vcf \
  --remove-indels \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode \
  --recode-INFO-all \
  --out data/pasteuria_population_filtered

# This creates: data/pasteuria_population_filtered.recode.vcf
# Rename it for clarity
mv data/pasteuria_population_filtered.recode.vcf data/pasteuria_population_filtered.vcf
```

**What these vcftools options do:**
- `--remove-indels`: Removes all insertion/deletion variants
- `--min-alleles 2 --max-alleles 2`: Keeps only biallelic sites (exactly 2 alleles)
- `--recode --recode-INFO-all`: Creates a new VCF file with all INFO fields preserved
- `--out`: Specifies the output file prefix

**Alternative: Using bcftools (more concise)**

```bash
# Filter with bcftools
bcftools view -v snps -m2 -M2 data/pasteuria_population.vcf -o data/pasteuria_population_filtered.vcf

# Options: -v snps (SNPs only), -m2 -M2 (min/max 2 alleles)
```

**Verify the filtering:**

```bash
# Check the original file
echo "Original file:"
vcftools --vcf data/pasteuria_population.vcf

# Check the filtered file
echo "Filtered file (SNPs only, biallelic):"
vcftools --vcf data/pasteuria_population_filtered.vcf
```

**This is your final filtered file for analysis:**
- Use `data/pasteuria_population_filtered.vcf` for all subsequent analyses

---

### 3.4 Extract Sample Names

It's useful to document which samples you're working with.

**Method 1: Using grep and text processing**

```bash
# Extract sample names from the VCF header (starting at column 10)
grep "^#CHROM" data/pasteuria_population_filtered.vcf | tr '\t' '\n' | tail -n +10 > results/sample_names.txt

# View the samples
cat results/sample_names.txt

# Count samples
wc -l results/sample_names.txt
```

**Alternative: Using bcftools**

```bash
# bcftools can query sample names directly
bcftools query -l data/pasteuria_population_filtered.vcf > results/sample_names.txt

# View and count
cat results/sample_names.txt
wc -l results/sample_names.txt
```

---

### 3.5 Summary of Data Exploration

Before moving on, record these key values (you'll need them later):

- **Number of samples (n)**: _______
- **Number of variant sites after filtering (S)**: _______
- **Total sequence length (L)**: _______ bp (ask your instructor)

**Note:** The sequence length represents the total genomic region analyzed, not just the number of variant sites.

---

## Section 4: Calculate Genetic Diversity Statistics

Now we'll calculate the two key diversity measures: S and π.

### 4.1 Calculate S (Number of Segregating Sites)

This is the simplest statistic - just count the polymorphic sites in our filtered VCF!

**Method 1: Using vcftools (reports automatically)**

```bash
# vcftools reports the number of sites automatically
vcftools --vcf data/pasteuria_population_filtered.vcf
```

This will print output like: "After filtering, kept 1234 out of a possible 1234 Sites"

**Method 2: Count lines directly with grep**

```bash
# Count lines in the VCF (excluding header lines that start with #)
grep -v "^#" data/pasteuria_population_filtered.vcf | wc -l

# Save S to results file
echo "S = $(grep -v '^#' data/pasteuria_population_filtered.vcf | wc -l)" > results/diversity_statistics.txt
```

**Alternative: Using bcftools**

```bash
# bcftools can also count variants
bcftools view -H data/pasteuria_population_filtered.vcf | wc -l
```

**Write down your value for S:** _______

**What does this tell us?**
- S is the total number of variable positions in our dataset
- Higher S = more mutations have occurred and been retained in the population
- S increases with mutation rate, population size, and sample size

---

### 4.2 Calculate π (Nucleotide Diversity) with VCFtools

VCFtools calculates π (nucleotide diversity) as the average proportion of nucleotide differences between all pairs of sequences.

**Calculate genome-wide π:**

```bash
# Calculate nucleotide diversity across the entire region
# Use a very large window size (larger than your total sequence length)
vcftools --vcf data/pasteuria_population_filtered.vcf \
  --window-pi 10000000 \
  --out results/pasteuria_pi

# View the result - the PI column (column 5) is your answer
cat results/pasteuria_pi.windowed.pi
```

The output has columns: CHROM, BIN_START, BIN_END, N_VARIANTS, PI
- **PI** (column 5) is your mean nucleotide diversity (π)

**Save π to your results:**

```bash
# Extract just the PI value and save it
tail -n 1 results/pasteuria_pi.windowed.pi | awk '{print "pi = " $5}' >> results/diversity_statistics.txt
```

**Record your π value:** _______

**What does this tell us?**
- π represents the average genetic distance between any two sequences
- π accounts for allele frequencies (common variants contribute more)
- Typical values: humans ~0.001, Drosophila ~0.01
- Higher π = more genetic diversity in the population

### 4.3 Visualize Diversity Across the Genome (Optional)

Instead of just calculating genome-wide averages, we can look at how diversity varies across the genome using sliding windows.

```bash
# Calculate pi in 10kb windows
vcftools --vcf data/pasteuria_population_filtered.vcf \
  --window-pi 10000 \
  --out results/pasteuria_pi_10kb

# View the first few windows
head -n 20 results/pasteuria_pi_10kb.windowed.pi
```

This shows you π calculated in sliding windows across the genome. You could plot this data to visualize:
- Regions of high vs. low diversity
- Potential selective sweeps (very low diversity regions)
- Highly variable regions

---

## Section 5: Interpret Results

### 5.1 Summary of Results

View all your calculated statistics:

```bash
# View the summary
cat results/diversity_statistics.txt
```

Create your own summary table:

| Statistic | Your Value | Units | What it means |
|-----------|------------|-------|---------------|
| S (segregating sites) | | sites | Number of polymorphic positions |
| π (pi) | | per site | Average pairwise differences |

---

## Section 6: Reflection and Discussion

### Discussion Questions

1. **Sample size effects:**
   - How might the number of samples (n) affect S?
   - Would π change if you sampled more individuals from the same population?
   - Why is it important to report sample size with diversity statistics?

2. **Biological interpretation:**
   - What might high vs. low genetic diversity tell us about this *Pasteuria* population?
   - How could you use these statistics to compare different populations or species?
   - What evolutionary processes increase diversity? Which ones decrease it?

3. **Methodological considerations:**
   - What assumptions does π make about the data?
   - What potential biases might exist in the VCF data (sequencing errors, missing data)?
   - Why did we filter out indels and multiallelic sites? Was it really necessary for the statistics we calculated?

---

## Section 7: Clean Up and Organize

Good practice is to document your analysis and organize your results:

```bash
# View all result files
ls -lh results/

# Remove any temporary files
rm -f scratch/*
```

**Your results directory should contain:**
- `diversity_statistics.txt` - Summary of S and π values
- `sample_names.txt` - List of samples analyzed
- `pasteuria_pi.windowed.pi` - Genome-wide π calculation
- `pasteuria_pi_10kb.windowed.pi` - Windowed π values (if you ran the optional section)

---

## Additional Resources

### Further Reading

- **VCF Format Specification**: https://samtools.github.io/hts-specs/VCFv4.2.pdf
- **VCFtools Documentation**: https://vcftools.github.io/
- **VCFtools Manual**: Run `man vcftools` in the terminal


### Commands Summary

```bash
# Quick reference for diversity calculations with vcftools

# 1. Filter VCF to biallelic SNPs only
vcftools --vcf data/file.vcf \
  --remove-indels \
  --min-alleles 2 --max-alleles 2 \
  --recode --recode-INFO-all \
  --out data/file_filtered

# 2. Count segregating sites
vcftools --vcf data/file.vcf
# Or: grep -v "^#" data/file.vcf | wc -l

# 3. Calculate pi genome-wide (assuming genome is less than 10 MB)
vcftools --vcf data/file.vcf --window-pi 10000000 --out results/pi

# 4. Calculate pi in sliding windows
vcftools --vcf data/file.vcf --window-pi 10000 --out results/pi_windows
```

---

## Assessment Questions

To check your understanding:

1. **What does S tell us that π doesn't, and vice versa?**
   - Think about rare vs. common variants

2. **If you sampled 5 individuals instead of 20, would you expect S to increase or decrease? What about π?**
   - Consider which statistic is more affected by sample size

3. **Why might different regions of the genome show different levels of diversity?**
   - Think about selection, recombination, mutation rate

4. **Looking at your results, would you say this *Pasteuria* population has high or low genetic diversity?**
   - What factors might influence diversity in a bacterial population?

---

**End of Practical Exercise**

**Remember to save all your results and notes! You may need them for your exam or future analyses.**

---

## Homework Assignment: Exploring Window Size Effects on Diversity Estimation

### Objective

Investigate how window size affects nucleotide diversity (π) estimates across the genome and determine if there's an optimal window size for analyzing this dataset.

### Task

**Part 1: Calculate π with Multiple Window Sizes**

Calculate nucleotide diversity across the genome using VCFtools with several different sliding window sizes. Choose at least 4-5 different window sizes that span a reasonable range (e.g., from small windows like 1kb to large windows like 100kb).

Consider:
- What range of window sizes makes sense given your genome/region size?
- What biological features might you want to capture with different window sizes?

**Part 2: Visualize the Results in R**

Load all of your windowed π results into R and create a single publication-quality plot with:
- **X-axis**: Genomic position (in bp or kb)
- **Y-axis**: Nucleotide diversity (π)
- **Multiple lines**: One line for each window size, each in a different color
- **Line styles**:
  - Your preferred/optimal window size should be plotted as a **solid line**
  - All other window sizes should be plotted as **dashed lines**
- **Legend**: Clearly showing which color corresponds to which window size
- **Axis labels**: Properly labeled with units
- **Title**: Descriptive title

### Deliverables

Submit the following two files:

1. **Your R script** - The complete R code used to:
   - Load the VCFtools output files
   - Create the plot
   - Export the PNG file

2. **Your plot** - A high-resolution PNG file showing:
   - All window sizes in different colors
   - Your preferred window size as a solid line
   - Other window sizes as dashed lines
   - A legend identifying window sizes by color

### Tips

- Start with a small range of window sizes to test your approach
- Think about how to handle overlapping windows (VCFtools has options for this)
- Consider using the window midpoint as the x-coordinate for plotting
- Make sure your colors are distinguishable and colorblind-friendly
- Test your R code on one file first before loading all of them
- Include comments in your R script explaining what each section does

**End of Practical Exercise**
