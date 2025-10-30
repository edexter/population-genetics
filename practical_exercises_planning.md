# Bioinformatics Practicals Plan
## Masters Course: Bioinformatics and Population Genetics

**Course Duration:** 4 weeks  
**Number of Practicals:** 4  
**Computational Environment:** Student laptops (local processing)  
**Time Constraint:** Each practical should require less than in 10-15 minutes of wall-time CPU processing

---

## Overall Learning Goals

By the end of these practicals, students will be able to:
1. Navigate the complete bioinformatics pipeline from raw data to publication-ready analyses
2. Work confidently with standard file formats (FASTQ, SAM/BAM, VCF)
3. Perform essential population genetics analyses
4. Interpret results in a biological context
5. Use command-line tools efficiently

---

## Dataset Strategy

### Primary Dataset: *E. coli* strains
- **Number of samples:** 10-12 samples from 3 distinct lineages/pathotypes (3-4 samples per lineage)
- **Reference genome:** *E. coli* K-12 MG1655 (~4.6 Mb)
- **Sequencing data:** Illumina paired-end reads, 150bp
- **Read depth:** ~50,000 read pairs per sample (~20-30x coverage)
- **Expected variants:** 5,000-15,000 SNPs depending on divergence

### Why *E. coli*?
- Small genome (fast processing on laptops)
- Well-annotated and medically relevant
- Clear population structure (different pathotypes)
- Abundant public data available
- Biologically meaningful for students

### Alternative Datasets to Consider:
- **Mycobacterium tuberculosis** outbreak simulation (smaller genome, public health angle)
- **Saccharomyces cerevisiae** chromosome (domestication biology, well-studied)
- **Pre-made RADseq dataset** (skip to population genetics if time is limited)

---

## Practical 1: Quality Control and Read Alignment

### Learning Objectives
- Download reference genomes from NCBI
- Assess raw sequencing read quality
- Align reads to a reference genome
- Understand SAM/BAM format and basic manipulation
- Interpret alignment statistics

### Software Required
- `wget` or `curl` (downloading data)
- FastQC (optional - for QC visualization)
- BWA (read alignment)
- samtools (BAM manipulation)

### Dataset for Practical 1
- Reference: *E. coli* K-12 genome (download from NCBI)
- Reads: 2-3 sample FASTQ files (paired-end)
- Total size: ~500 MB compressed

### Practical Steps
1. **Download reference genome** from NCBI using command line
   - Learn to navigate NCBI FTP site
   - Use `wget` or `curl`
   - Verify file integrity

2. **Index reference genome** with BWA
   ```bash
   bwa index reference.fasta
   ```

3. **(Optional) Run FastQC** on one sample
   - Visualize read quality
   - Identify any quality issues
   - Discuss what "good" data looks like

4. **Align reads with BWA-MEM**
   ```bash
   bwa mem reference.fasta reads_R1.fastq reads_R2.fastq > alignment.sam
   ```

5. **Convert SAM to BAM, sort, and index**
   ```bash
   samtools view -b alignment.sam > alignment.bam
   samtools sort alignment.bam -o alignment_sorted.bam
   samtools index alignment_sorted.bam
   ```

6. **View alignment statistics**
   ```bash
   samtools flagstat alignment_sorted.bam
   samtools depth alignment_sorted.bam | head
   ```

### Expected Computation Time
- Download: 2-3 minutes
- Indexing: <1 minute
- Alignment (per sample): 3-5 minutes
- SAM/BAM conversion: 1-2 minutes
- **Total: ~10 minutes per sample**

### Deliverables
- Sorted, indexed BAM files for each sample
- Alignment statistics summary
- Understanding of the alignment pipeline

### Key Concepts to Emphasize
- File formats: FASTQ, SAM, BAM
- What alignment means biologically
- Coverage and depth
- Why we sort and index BAM files

---

## Practical 2: Variant Calling and VCF Manipulation

### Learning Objectives
- Call variants from aligned reads
- Understand VCF format structure
- Filter variants by quality metrics
- Extract and manipulate variant data
- Calculate basic variant statistics

### Software Required
- bcftools (variant calling and manipulation)
- VCFtools (filtering and statistics)
- Basic Unix tools (grep, awk, wc)

### Dataset for Practical 2
- Input: Sorted BAM files from Practical 1 (3-5 samples)
- Reference genome (same as Practical 1)

### Practical Steps

1. **Call variants with bcftools**
   ```bash
   bcftools mpileup -f reference.fasta *.bam | bcftools call -mv -Oz -o variants_raw.vcf.gz
   ```

2. **Examine the VCF format**
   ```bash
   zcat variants_raw.vcf.gz | head -50
   zcat variants_raw.vcf.gz | grep -v "^#" | head -10
   ```
   - Discuss header lines
   - Explain each column
   - Show how genotypes are encoded

3. **Filter variants by quality**
   ```bash
   bcftools view -i 'QUAL>20 && DP>10' variants_raw.vcf.gz -Oz -o variants_filtered.vcf.gz
   ```

4. **Calculate basic statistics**
   ```bash
   bcftools stats variants_filtered.vcf.gz > variant_stats.txt
   ```
   - Count total SNPs
   - Calculate transition/transversion ratio
   - Examine variant density

5. **Extract specific information**
   ```bash
   # Count SNPs per sample
   bcftools query -l variants_filtered.vcf.gz
   
   # Extract SNP positions
   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' variants_filtered.vcf.gz > snp_list.txt
   ```

### Expected Computation Time
- Variant calling: 3-5 minutes
- Filtering and statistics: 1-2 minutes
- **Total: ~5-7 minutes**

### Deliverables
- Filtered VCF file with high-quality SNPs
- Summary statistics document
- Understanding of variant filtering rationale

### Key Concepts to Emphasize
- What is a variant call?
- Quality metrics (QUAL, depth, mapping quality)
- Why filtering is necessary
- VCF format structure
- Difference between SNPs and indels

---

## Practical 3: Population Genetics Analysis

### Learning Objectives
- Convert between genomic file formats
- Calculate population genetic statistics (π, Fst)
- Perform Principal Component Analysis (PCA)
- Infer population structure with ADMIXTURE
- Interpret population genetics results

### Software Required
- PLINK (file conversion, PCA)
- VCFtools (diversity statistics)
- ADMIXTURE (population structure)
- R (visualization - provide scripts)

### Dataset for Practical 3
**Option A:** Use VCF from Practical 2 (if sufficient samples)  
**Option B:** Provide pre-made VCF with 10-12 samples from 3 populations (recommended)

### Practical Steps

1. **Convert VCF to PLINK format**
   ```bash
   plink --vcf variants_filtered.vcf.gz --make-bed --out variants --allow-extra-chr
   ```

2. **Calculate nucleotide diversity (π)**
   ```bash
   vcftools --gzvcf variants_filtered.vcf.gz --window-pi 10000 --out diversity
   ```

3. **Calculate pairwise Fst between populations**
   - Requires population assignment file
   ```bash
   vcftools --gzvcf variants_filtered.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out fst_pop1_pop2
   ```

4. **Run PCA**
   ```bash
   plink --bfile variants --pca --out pca_results
   ```

5. **Run ADMIXTURE for K=2 and K=3**
   ```bash
   admixture variants.bed 2
   admixture variants.bed 3
   ```

6. **Visualize results in R** (provide pre-made scripts)
   - PCA plot
   - ADMIXTURE barplot
   - Diversity across genome

### Expected Computation Time
- File conversion: <1 minute
- Statistics calculation: 1-2 minutes
- PCA: <1 minute
- ADMIXTURE: 2-3 minutes per K
- Visualization: 1-2 minutes
- **Total: ~8-10 minutes**

### Deliverables
- Summary statistics table (π, Fst)
- PCA plot showing population structure
- ADMIXTURE barplot
- Written interpretation of results

### Key Concepts to Emphasize
- What nucleotide diversity tells us
- Interpreting Fst values (population differentiation)
- PCA as dimensionality reduction
- What ADMIXTURE ancestry coefficients mean
- Relating statistics to biology (selection, drift, gene flow)

### Important Note
For this practical to work well, you need samples from distinct populations. If using student-generated data from Practicals 1-2, you may need to provide a pre-made dataset with clear population structure.

---

## Practical 4: Phylogenetic Analysis

### Learning Objectives
- Prepare sequence data for phylogenetic analysis
- Build phylogenetic trees using maximum likelihood
- Visualize and annotate phylogenetic trees
- Interpret evolutionary relationships
- Relate phylogenies to population structure

### Software Required
**Option A (SNP-based - recommended):**
- vcf2phylip or custom script
- IQ-TREE or FastTree
- FigTree or R with ggtree

**Option B (Gene-based - simpler alternative):**
- MUSCLE or MAFFT (alignment)
- FastTree (fast ML tree)
- FigTree (visualization)

### Dataset for Practical 4

**Option A (SNP-based):** Use VCF from Practicals 2-3  
**Option B (Gene-based):** 10-15 bacterial 16S rRNA sequences from NCBI

### Practical Steps (Option A - SNP-based)

1. **Convert VCF to phylip alignment format**
   ```bash
   python vcf2phylip.py -i variants_filtered.vcf.gz
   ```

2. **Build maximum likelihood tree**
   ```bash
   iqtree -s variants.phy -m GTR+ASC -bb 1000 -nt AUTO
   ```
   OR (faster)
   ```bash
   FastTree -nt -gtr variants.phy > tree.newick
   ```

3. **Visualize tree in FigTree**
   - Open tree file
   - Root tree (if appropriate)
   - Add sample labels
   - Color by population
   - Export figure

4. **Interpret tree topology**
   - Which samples are most closely related?
   - Does the tree match population assignments?
   - Are there any surprises?

### Practical Steps (Option B - Gene-based)

1. **Download 16S rRNA sequences from NCBI**
   - Provide accession numbers or pre-download

2. **Align sequences with MUSCLE**
   ```bash
   muscle -in sequences.fasta -out aligned.fasta
   ```

3. **Build tree with FastTree**
   ```bash
   FastTree -nt aligned.fasta > tree.newick
   ```

4. **Visualize and interpret** (same as above)

### Expected Computation Time
- Alignment/conversion: 1-2 minutes
- Tree building: 3-5 minutes (faster with FastTree)
- Visualization: 2-3 minutes
- **Total: ~8-10 minutes**

### Deliverables
- Phylogenetic tree figure (publication quality)
- Written interpretation of evolutionary relationships
- Comparison to population structure from Practical 3

### Key Concepts to Emphasize
- What a phylogenetic tree represents
- Maximum likelihood vs. other methods
- Bootstrap support values
- Rooting trees
- Gene trees vs. species trees
- Connection between phylogeny and population genetics

---

## Software Installation Plan

### Before Practical 1
Students should install:
- BWA
- samtools
- FastQC (optional)

### Before Practical 2
Students should install:
- bcftools
- VCFtools

### Before Practical 3
Students should install:
- PLINK
- ADMIXTURE
- R (with basic packages)

### Before Practical 4
Students should install:
- IQ-TREE or FastTree
- FigTree
- (Optional) MUSCLE/MAFFT

### Installation Resources
- Provide conda environment YAML file for easy installation
- Test all installations on different OS (Mac, Linux, Windows with WSL)
- Have backup: cloud instances or JSLinux if local installation fails

---

## Data Preparation Checklist

### For Instructor (You)

**Before Course Starts:**
- [ ] Download or simulate *E. coli* whole genome sequencing data
- [ ] Downsample reads to appropriate coverage (~50k read pairs per sample)
- [ ] Test complete pipeline on your machine
- [ ] Time each practical step
- [ ] Create population assignment files for Practical 3
- [ ] Prepare R visualization scripts
- [ ] Write detailed protocols for each practical
- [ ] Create answer keys

**Datasets to Provide:**
- [ ] Reference genome FASTA
- [ ] 10-12 FASTQ file pairs (compressed)
- [ ] (Optional) Pre-made VCF for Practical 3 if timeline is tight
- [ ] (Optional) 16S sequences for Practical 4 Option B
- [ ] Population metadata file
- [ ] Sample information sheet

**Total Data Size:** ~2-3 GB compressed

---

## Contingency Plans

### If Computation is Too Slow
- Skip Practical 1, provide pre-made BAM files
- Skip Practical 2, provide pre-made VCF
- Focus hands-on time on Practicals 3-4 (population genetics core)

### If Installation Issues Arise
- Use conda environments (easiest)
- Have cloud computing backup (AWS, Google Cloud)
- Use browser-based tools where available
- Pair students (one runs while both observe)

### If Dataset is Too Large
- Use single chromosome instead of whole genome
- Use RADseq-style data (fewer loci)
- Further downsample reads
- Use bacterial genome instead of eukaryote

---

## Assessment Ideas

### Formative Assessment (During Practicals)
- Checkpoint questions at each step
- Predict outcomes before running commands
- Interpret output files
- Troubleshooting exercises

### Summative Assessment (End of Practicals)
- Written report interpreting their results
- Present findings to class (5 min presentation)
- Answer biological questions based on their data
- Design their own analysis pipeline for a new question

### Example Questions
1. "Based on your PCA and ADMIXTURE results, how many distinct populations are in your dataset? What biological factors might explain this structure?"

2. "Your Fst value between populations A and B is 0.15. What does this tell you about gene flow?"

3. "How does the phylogenetic tree from Practical 4 compare to the population structure from Practical 3? Are they consistent?"

---

## Timeline Recommendation

### Week 1: Linux Basics + Practical 1
- Days 1-2: Linux command line training
- Days 3-4: Practical 1 (QC and Alignment)

### Week 2: Practical 2 + Population Genetics Theory
- Days 1-2: Practical 2 (Variant Calling)
- Days 3-4: Theory lectures on population genetics

### Week 3: Practical 3
- Days 1-3: Practical 3 (Population Genetics Analysis)
- Day 4: Interpretation and discussion

### Week 4: Practical 4 + Synthesis
- Days 1-2: Practical 4 (Phylogenetics)
- Days 3-4: Student presentations and wrap-up

---

## Resources and References

### Software Documentation
- BWA: http://bio-bwa.sourceforge.net/
- samtools: http://www.htslib.org/
- bcftools: http://samtools.github.io/bcftools/
- PLINK: https://www.cog-genomics.org/plink/
- ADMIXTURE: https://dalexander.github.io/admixture/
- IQ-TREE: http://www.iqtree.org/

### Tutorials and Guides
- GATK Best Practices (for reference): https://gatk.broadinstitute.org/
- Heng Li's blog (samtools author): https://lh3.github.io/
- Population genetics in R: https://grunwaldlab.github.io/Population_Genetics_in_R/

### Public Datasets
- NCBI SRA: https://www.ncbi.nlm.nih.gov/sra
- EBI ENA: https://www.ebi.ac.uk/ena
- *E. coli* reference genomes: https://www.ncbi.nlm.nih.gov/genome/167

---

## Notes and Considerations

### Pedagogical Approach
- Start with big picture (what are we trying to learn?)
- Show example outputs BEFORE students generate them
- Emphasize biological interpretation over technical details
- Connect each practical to real research questions
- Encourage troubleshooting and problem-solving

### Common Student Struggles
- File paths (absolute vs. relative) - address early
- File format confusion - show visual examples
- Long command lines - provide templates
- Interpreting statistics - use analogies and visualizations
- Patience during computation - have discussion questions ready

### Future Improvements
- Add practical on RNA-seq analysis
- Include more visualization in R
- Add GWAS practical if time permits
- Consider adding ancient DNA or metagenomics module

---

## Success Criteria

Students should be able to:
- [ ] Execute a complete bioinformatics pipeline independently
- [ ] Troubleshoot common command-line errors
- [ ] Interpret population genetics statistics correctly
- [ ] Make publication-quality figures
- [ ] Connect analyses to biological questions
- [ ] Design appropriate analyses for new research questions

---

**Last Updated:** [Date]  
**Instructor:** [Your Name]  
**Course:** Bioinformatics and Population Genetics  
**Institution:** [University Name]