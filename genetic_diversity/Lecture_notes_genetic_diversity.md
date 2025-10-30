# Measuring Genetic Diversity in Population Genomics
## Lecture Notes

---

## Introduction

### Why Measure Genetic Diversity?

Genetic diversity is fundamental to understanding:
- **Evolution**: How populations adapt and change over time
- **Conservation**: Identifying populations at risk due to low diversity
- **Health**: Understanding disease susceptibility and resistance
- **Breeding**: Selecting for desired traits in agriculture

### What We'll Learn Today

In this lecture, we will cover four main statistics used to quantify genetic diversity:
1. Number of Segregating Sites (S)
2. Nucleotide Diversity (π)
3. Heterozygosity (He)
4. Watterson's Theta (θw)

We will also explore how these statistics relate to each other and what they tell us about population history.

### Our Example Dataset

Throughout this lecture, we'll work with a simple example:
- **5 DNA sequences** (samples from a population)
- **10 base pairs** long (aligned)
- We'll calculate each statistic by hand

**Example Alignment:**
```
Seq1: A T G C A T G C A T
Seq2: A T G C G T G C A T
Seq3: A C G C A T G C A T
Seq4: A T G C A T A C A T
Seq5: A T G T A T G C A T
```

---

## 1. Number of Segregating Sites (S)

### Definition

The **number of segregating sites** (S) is simply the count of polymorphic (variable) positions in an alignment. A site is polymorphic if it has more than one nucleotide variant in the sample.

### Hand Calculation

Looking at our example alignment:

```
Position: 1 2 3 4 5 6 7 8 9 10
Seq1:     A T G C A T G C A T
Seq2:     A T G C G T G C A T
Seq3:     A C G C A T G C A T
Seq4:     A T G C A T A C A T
Seq5:     A T G T A T G C A T
```

**Step 1:** Examine each position
- Position 1: A, A, A, A, A → **Not variable** (monomorphic)
- Position 2: T, T, C, T, T → **Variable** (T and C)
- Position 3: G, G, G, G, G → **Not variable**
- Position 4: C, C, C, C, T → **Variable** (C and T)
- Position 5: A, G, A, A, A → **Variable** (A and G)
- Position 6: T, T, T, T, T → **Not variable**
- Position 7: G, G, G, A, G → **Variable** (G and A)
- Position 8: C, C, C, C, C → **Not variable**
- Position 9: A, A, A, A, A → **Not variable**
- Position 10: T, T, T, T, T → **Not variable**

**Step 2:** Count variable sites

S = 4 (positions 2, 4, 5, and 7 are variable)

### Interpretation

- S is the simplest measure of genetic variation
- **Important caveat**: S is dependent on sample size—larger samples will detect more rare variants
- S increases with:
  - Higher mutation rate
  - Larger effective population size
  - Larger sample size

### Typical Values

- Highly variable depending on organism, genomic region, and sample size
- Low diversity regions: S might be 0 or very small
- High diversity regions: S can be quite large relative to sequence length

---

## 2. Nucleotide Diversity (π)

### Definition

**Nucleotide diversity** (π, pi) is the average number of nucleotide differences per site between two randomly chosen sequences from the population. It's a measure of genetic variation that accounts for allele frequencies.

### Formula

$$\pi = \frac{\sum \text{(pairwise differences)}}{\text{number of comparisons} \times \text{sequence length}}$$

Alternatively:

$$\pi = \frac{1}{\binom{n}{2} L} \sum_{ {i < j} } d_{ij}$$

Where:
- n = number of sequences
- L = sequence length
- dij = number of differences between sequences i and j

### Hand Calculation

Using our example alignment with 5 sequences and 10 sites:

**Step 1:** Calculate all pairwise comparisons

For 5 sequences, we have $\binom{5}{2} = 10$ pairs:

| Pair | Seq1 vs Seq2 | Seq1 vs Seq3 | Seq1 vs Seq4 | Seq1 vs Seq5 | Seq2 vs Seq3 |
|------|--------------|--------------|--------------|--------------|--------------|
| Differences | 1 | 1 | 1 | 1 | 2 |

| Pair | Seq2 vs Seq4 | Seq2 vs Seq5 | Seq3 vs Seq4 | Seq3 vs Seq5 | Seq4 vs Seq5 |
|------|--------------|--------------|--------------|--------------|--------------|
| Differences | 2 | 2 | 2 | 2 | 2 |

**Step 2:** Count differences for each pair

Let's examine Seq1 vs Seq2:
```
Seq1: A T G C A T G C A T
Seq2: A T G C G T G C A T
Diff:         *           
```
Differences = 1 (at position 5)

Continuing for all pairs:
- Seq1 vs Seq2: 1 difference
- Seq1 vs Seq3: 1 difference (position 2)
- Seq1 vs Seq4: 1 difference (position 7)
- Seq1 vs Seq5: 1 difference (position 4)
- Seq2 vs Seq3: 2 differences (positions 2, 5)
- Seq2 vs Seq4: 2 differences (positions 5, 7)
- Seq2 vs Seq5: 2 differences (positions 4, 5)
- Seq3 vs Seq4: 2 differences (positions 2, 7)
- Seq3 vs Seq5: 2 differences (positions 2, 4)
- Seq4 vs Seq5: 2 differences (positions 4, 7)

**Step 3:** Sum and divide

Total differences = 1 + 1 + 1 + 1 + 2 + 2 + 2 + 2 + 2 + 2 = 16

$$\pi = \frac{16}{10 \text{ pairs} \times 10 \text{ sites}} = \frac{16}{100} = 0.16$$

### Interpretation

- π represents the **per-site diversity**
- π = 0.16 means that, on average, 16% of nucleotide positions differ between two random sequences
- Under neutrality, π is an estimate of θ = 4Neμ, where:
  - Ne = effective population size
  - μ = mutation rate per site per generation
- π weights by allele frequency: common variants contribute more than rare variants

### Typical Values

- In humans: π ≈ 0.001 (0.1%) across the genome
- In Drosophila: π ≈ 0.01 (1%)
- Highly variable depending on:
  - Species effective population size
  - Mutation rate
  - Genomic region (coding vs non-coding)
  - Population history (bottlenecks reduce π)

---

## 3. Heterozygosity (He)

### Definition

**Expected heterozygosity** (He), also called **gene diversity**, is the probability that two randomly chosen alleles at a locus are different. It's calculated based on allele frequencies and represents the expected proportion of heterozygotes under Hardy-Weinberg equilibrium.

### Formula

For a single site:

$$H_e = 1 - \sum_{i} p_i^2$$

Where pi is the frequency of allele i at that site.

For multiple sites, average He across all sites (or across variable sites only).

### Hand Calculation

Using our example alignment:

**Step 1:** Calculate allele frequencies at each variable site

Position 2:
- T: 4 out of 5 (frequency = 0.8)
- C: 1 out of 5 (frequency = 0.2)

$$H_{e,2} = 1 - (0.8^2 + 0.2^2) = 1 - (0.64 + 0.04) = 1 - 0.68 = 0.32$$

Position 4:
- C: 4 out of 5 (frequency = 0.8)
- T: 1 out of 5 (frequency = 0.2)

$$H_{e,4} = 1 - (0.8^2 + 0.2^2) = 0.32$$

Position 5:
- A: 4 out of 5 (frequency = 0.8)
- G: 1 out of 5 (frequency = 0.2)

$$H_{e,5} = 1 - (0.8^2 + 0.2^2) = 0.32$$

Position 7:
- G: 4 out of 5 (frequency = 0.8)
- A: 1 out of 5 (frequency = 0.2)

$$H_{e,7} = 1 - (0.8^2 + 0.2^2) = 0.32$$

**Step 2:** Average across variable sites

$$\bar{H}_e = \frac{0.32 + 0.32 + 0.32 + 0.32}{4} = 0.32$$

Or, if averaging across all sites (including monomorphic):

$$\bar{H}_e = \frac{0.32 + 0.32 + 0.32 + 0.32}{10} = 0.128$$

### Interpretation

- He measures the **expected genetic diversity at a locus**
- Higher He indicates more balanced allele frequencies
- He = 0: No variation (all individuals have the same allele)
- He = 0.5: Maximum diversity for a biallelic locus (both alleles at 50%)
- Related to π, but calculated per-site then averaged
- Both weight by allele frequency

### Note on Observed Heterozygosity (Ho)

If working with diploid genotype data, you can also calculate:
- **Ho** (observed heterozygosity): The actual proportion of heterozygotes in the sample
- Comparing Ho and He can reveal inbreeding or population structure

---

## 4. Watterson's Theta (θw)

### Definition

**Watterson's theta** (θw) is an estimate of the population mutation parameter θ = 4Neμ based on the number of segregating sites (S). Unlike π, it treats all mutations equally regardless of their frequency.

### Formula

$$\theta_w = \frac{S}{a_n}$$

Where:
- S = number of segregating sites
- an = harmonic number = $\sum_{i=1}^{n-1} \frac{1}{i}$
- n = number of sequences (sample size)

### Why We Need It

The harmonic number correction (an) accounts for the fact that larger sample sizes will detect more rare variants. This makes θw comparable across studies with different sample sizes.

### Hand Calculation

Using our example with n = 5 sequences and S = 4:

**Step 1:** Calculate the harmonic number

$$a_n = \sum_{i=1}^{n-1} \frac{1}{i} = \sum_{i=1}^{4} \frac{1}{i}$$

$$a_5 = \frac{1}{1} + \frac{1}{2} + \frac{1}{3} + \frac{1}{4}$$

$$a_5 = 1 + 0.5 + 0.333 + 0.25 = 2.083$$

**Step 2:** Calculate θw

$$\theta_w = \frac{S}{a_n} = \frac{4}{2.083} = 1.92$$

**Step 3:** Convert to per-site estimate (optional)

To make it comparable to π, divide by sequence length:

$$\theta_w \text{ per site} = \frac{1.92}{10} = 0.192$$

### Interpretation

- θw is another estimator of θ = 4Neμ
- Under neutrality and constant population size, θw and π should be similar
- θw is based only on the **number** of segregating sites, not their frequencies
- More sensitive to rare variants than π
- Used in neutrality tests (e.g., Tajima's D)

### Typical Values

- Scale similar to π
- Can be higher or lower than π depending on allele frequency spectrum

---

## 5. Relationships Between Statistics

### Comparing π and θw

Both π and θw estimate the same underlying parameter (θ = 4Neμ) but in different ways:

| Statistic | Based on | Weighting | Sensitivity |
|-----------|----------|-----------|-------------|
| π | Pairwise differences | Frequency-weighted | Common variants |
| θw | Number of segregating sites | Unweighted | Rare variants |

### What Does the Difference Mean?

The relationship between π and θw can reveal information about population history and selection:

**1. π ≈ θw**
- Consistent with neutrality and constant population size
- Allele frequency spectrum is as expected under neutral evolution

**2. π > θw**
- Excess of **intermediate-frequency variants**
- Possible explanations:
  - Balancing selection
  - Population bottleneck (recent)
  - Population structure

**3. π < θw**
- Excess of **rare variants**
- Possible explanations:
  - Recent population expansion
  - Purifying selection
  - Selective sweep (directional selection)

### Tajima's D (Conceptual Overview)

**Tajima's D** is a statistical test that formalizes the comparison between π and θw:

$$D = \frac{\pi - \theta_w}{\text{standard error}}$$

- **D = 0**: Consistent with neutral evolution and constant population size
- **D > 0**: π > θw (excess intermediate-frequency variants)
- **D < 0**: π < θw (excess rare variants)

**Important:** We won't calculate Tajima's D by hand (requires complex variance calculations), but understanding what the difference between π and θw means is crucial for interpreting diversity patterns.

### Connection to Heterozygosity

He and π are closely related:
- Both are frequency-weighted measures
- He is calculated per-site then averaged
- π is calculated from all pairwise comparisons
- In practice, they often give similar results
- He is more commonly used when thinking about single loci
- π is more commonly used for sequence-level diversity

### Our Example Results

From our calculations:
- S = 4
- π = 0.16
- He = 0.32 (per variable site) or 0.128 (per all sites)
- θw = 0.192 (per site)

In our example, π and θw are similar (0.16 vs 0.192), suggesting a relatively neutral allele frequency spectrum. The small difference could be due to sampling variation with such a small dataset.

---

## Summary

### Key Takeaways

1. **Different statistics capture different aspects of genetic diversity:**
   - **S**: Raw count of variable sites (sample-size dependent)
   - **π**: Average pairwise differences (frequency-weighted)
   - **He**: Expected heterozygosity (per-site diversity)
   - **θw**: Diversity based on segregating sites (corrected for sample size)

2. **All statistics increase with:**
   - Higher mutation rate (μ)
   - Larger effective population size (Ne)
   - For some statistics: larger sample size

3. **Comparing statistics reveals population history:**
   - π vs θw can indicate selection or demographic changes
   - Understanding what causes differences is key to interpretation

4. **Sample size matters:**
   - Always report sample size with diversity statistics
   - Use θw when comparing across different sample sizes
   - Larger samples detect more rare variants

### When to Use Which Statistic?

- **S**: Quick overview of variation, but not standardized
- **π**: Most commonly reported sequence diversity measure
- **He**: Comparing diversity at specific loci, especially with genotype data
- **θw**: When sample sizes differ across studies, or for neutrality tests

### Next Steps: The Practical

In the practical session, you will:
- Use real genomic data (VCF files or aligned sequences)
- Calculate these statistics using bioinformatic tools (e.g., vcftools, pixy)
- Interpret diversity patterns across a genome
- Compare diversity between populations or genomic regions

---

## Formulas Summary Sheet

### Number of Segregating Sites
$$S = \text{count of polymorphic sites}$$

### Nucleotide Diversity
$$\pi = \frac{\sum \text{pairwise differences}}{\binom{n}{2} \times L}$$

### Expected Heterozygosity
$$H_e = 1 - \sum_{i} p_i^2$$

### Watterson's Theta
$$\theta_w = \frac{S}{a_n}$$

where $a_n = \sum_{i=1}^{n-1} \frac{1}{i}$

### Tajima's D (concept only)
$$D = \frac{\pi - \theta_w}{\text{standard error}}$$

---

## Additional Notes

### Effective Population Size (Ne)

- Ne is the number of individuals that would give the same genetic diversity if the population were ideal
- Always smaller than census population size
- Affected by:
  - Unequal sex ratios
  - Variation in reproductive success
  - Population fluctuations over time

### Mutation Rate (μ)

- Varies across species and genomic regions
- Typical values: 10⁻⁸ to 10⁻⁹ per base pair per generation
- Higher in some organisms (e.g., RNA viruses)
- Can be estimated from:
  - Pedigree studies
  - Phylogenetic comparisons
  - Mutation accumulation experiments

### Hardy-Weinberg Equilibrium

Expected heterozygosity assumes Hardy-Weinberg equilibrium:
- Random mating
- No selection
- No migration
- No mutation (within a generation)
- Large population size

Deviations from these assumptions affect observed vs expected heterozygosity.

---

**End of Lecture Notes**

### This information needs to be added:
# Understanding θ = 4Neμ

## What is θ (theta)?

**θ (theta)** is the **population mutation parameter** - a fundamental quantity in population genetics that describes the **expected level of genetic diversity** in a population at equilibrium.

It combines two key biological parameters:
- **Ne** = effective population size
- **μ** = mutation rate per site per generation

---

## The Formula: θ = 4Neμ

### Breaking Down Each Component

**Ne (effective population size):**
- The number of individuals that **effectively contribute** to reproduction
- Always ≤ census population size
- Larger Ne → more genetic diversity

**μ (mutation rate):**
- Probability that a given nucleotide mutates per generation
- Typical values: 10⁻⁸ to 10⁻⁹ per base pair per generation
- Higher μ → more genetic diversity

**4:**
- Comes from diploid organisms (2 copies per individual × 2 for mutation-drift balance)
- For haploids, it would be θ = 2Neμ

---

## What Does θ Represent?

θ is the **expected genetic diversity under neutral evolution at mutation-drift equilibrium**.

### Mutation-Drift Balance

In a population at equilibrium:
- **Mutations** create new variation (adds diversity)
- **Genetic drift** removes variation through random sampling (reduces diversity)

θ represents the **balance point** where these two forces are equal.

---

## Why Do We Care About θ?

### 1. It's a Property of the Population (Not the Sample)

θ is the "true" diversity of the population. When we calculate π or θw from a sample, we're trying to **estimate** this underlying parameter.

Think of it like:
- **θ** = true population mean height
- **π or θw** = estimate of mean height from a sample

---

### 2. It Tells Us About Population History and Biology

If we estimate θ = 0.01 from our data, this tells us about:

**Effective population size:**
- If we know μ = 2.5 × 10⁻⁸, then:
$$N_e = \frac{\theta}{4\mu} = \frac{0.01}{4 \times 2.5 \times 10^{-8}} = 100,000$$

**Population history:**
- Large θ → large Ne or high μ → high genetic diversity
- Small θ → small Ne or low μ → low genetic diversity
- Changes in θ over time indicate population size changes

---

## Why Two Different Estimators (π and θw)?

Both π and θw are trying to estimate the **same thing** (θ), but they use different information from the data:

### π (Nucleotide Diversity)
**What it uses:** The **frequency** of variants (pairwise differences)

**How it estimates θ:**
- Counts all differences between sequences
- Weights by allele frequency
- Common variants contribute more

**When it equals θ:** Under neutral evolution at equilibrium

---

### θw (Watterson's Theta)
**What it uses:** The **number** of variants (segregating sites)

**How it estimates θ:**
- Counts how many variable sites exist
- Doesn't weight by frequency
- All variants contribute equally

**When it equals θ:** Under neutral evolution at equilibrium

---

## An Analogy

Imagine you want to know the **true average income (θ)** in a city:

**Method 1 (like π):**
- Survey people and calculate the mean of their actual incomes
- High earners pull the average up
- Sensitive to the distribution of incomes

**Method 2 (like θw):**
- Count how many different income brackets exist
- Use that count to estimate average income
- Less sensitive to how many people are in each bracket

Both methods estimate the same thing (average income), but use different information. If the income distribution is "normal" (neutral evolution), both give similar answers. If it's skewed, they differ.

---

## What Does Mutation-Drift Equilibrium Mean?

### At Equilibrium:
- Population size has been constant for many generations
- No selection acting on the sites
- Mutation introduces new variants at rate μ
- Drift removes variants at rate dependent on Ne
- These forces balance out

### The Math:
$$\text{Expected diversity} = \theta = 4N_e\mu$$

This is what π and θw are trying to estimate!

---

## When π ≈ θw ≈ θ

When a population is:
- At **mutation-drift equilibrium**
- Under **neutral evolution** (no selection)
- At **constant population size**

Then:
$$\pi \approx \theta_w \approx \theta = 4N_e\mu$$

All three are estimating the same parameter and should be similar.

---

## When π ≠ θw (Both Still Estimate θ, But Differently)

When they differ, it tells us something is violating our assumptions:

### Case 1: π > θw
**What happened:** Excess of intermediate-frequency alleles

**Possible causes:**
- Balancing selection (maintains multiple alleles)
- Population bottleneck (recent)
- Population structure

**Why they differ:**
- π weights by frequency → intermediate frequencies increase π
- θw just counts sites → doesn't care about frequency

Both are still estimating θ, but they're capturing different aspects of non-equilibrium dynamics.

---

### Case 2: π < θw
**What happened:** Excess of rare alleles

**Possible causes:**
- Population expansion (recent)
- Selective sweep (removed variation)
- Purifying selection (removes deleterious variants)

**Why they differ:**
- π weights by frequency → rare variants contribute little to π
- θw counts sites → rare variants count just as much

---

## Practical Example

### Population A:
- Ne = 10,000
- μ = 2.5 × 10⁻⁸ per site per generation

**Expected θ:**
$$\theta = 4 \times 10,000 \times 2.5 \times 10^{-8} = 0.001$$

If we sequence this population:
- π ≈ 0.001 (if at equilibrium)
- θw ≈ 0.001 (if at equilibrium)

---

### Population B: 
Same as A, but recently expanded from Ne = 1,000 to 10,000

**Current θ** (based on current Ne):
$$\theta = 4 \times 10,000 \times 2.5 \times 10^{-8} = 0.001$$

But our sample shows:
- π ≈ 0.0005 (lots of rare variants from expansion)
- θw ≈ 0.0012 (counts all those rare variants)
- π < θw → signature of population expansion!

Both are still trying to estimate θ, but they differ because the population hasn't reached equilibrium yet.

---

## Summary

**θ = 4Neμ is:**
- The **true, underlying genetic diversity** of the population
- Expected diversity at **mutation-drift equilibrium**
- A function of **population size** and **mutation rate**
- What we're trying to **estimate** from our sample data

**π and θw are:**
- Two different **estimators** of θ
- Use different information (frequencies vs counts)
- Should be equal at equilibrium
- Their difference reveals **evolutionary history**

**The key insight:**
When someone says "both π and θw estimate θ," they mean:
- There's a true population parameter (θ) we want to know
- We have two different ways to estimate it from sample data
- Under ideal conditions, both give us the same answer
- When they differ, we learn about what evolutionary forces are acting

It's like having two different thermometers (π and θw) measuring the same temperature (θ) - if they both read the same, conditions are normal. If they disagree, something interesting is happening!
