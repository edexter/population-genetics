# The Normal Distribution — A Primer for Population Genetics

## 1. Concept Overview

The normal distribution (also called the Gaussian distribution) models continuous variables that cluster symmetrically around a mean value. It's the most important distribution in statistics because many natural phenomena approximate it, and sums of random variables tend toward normality.

**Key characteristics:**
- Bell-shaped and symmetric around the mean
- Completely described by two parameters: mean (μ) and variance (σ²)
- Continuous (not discrete like binomial or Poisson)

**Population genetics examples:**
- Distribution of quantitative traits (height, weight) in a population
- Distribution of allele frequencies across many loci under drift
- Sampling distribution of genetic statistics (FST, heterozygosity estimates)
- Body size or other polygenic traits

## 2. The Formula

The probability density for a value *x* is:

$$f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}$$

**Where:**
* *x* = the value we're evaluating
* μ (mu) = the mean (center of the distribution)
* σ (sigma) = the standard deviation (spread)
* σ² = the variance
* *e* ≈ 2.71828 and π ≈ 3.14159

**Note:** For continuous distributions, we calculate probabilities over intervals, not at exact points.

## 3. The Standard Normal Distribution

The **standard normal** has μ = 0 and σ = 1, denoted as *Z* ~ N(0, 1).

Any normal variable *X* ~ N(μ, σ²) can be standardized:

$$Z = \frac{X - \mu}{\sigma}$$

This *Z*-score tells us how many standard deviations *X* is from the mean.

**Example:** If mean height is μ = 170 cm with σ = 10 cm, a person who is 185 cm tall has:

$$Z = \frac{185 - 170}{10} = 1.5$$

They are 1.5 standard deviations above average.

## 4. The Empirical Rule (68-95-99.7 Rule)

For any normal distribution:
- **68%** of values fall within ±1σ of the mean
- **95%** of values fall within ±2σ of the mean
- **99.7%** of values fall within ±3σ of the mean

This rule is invaluable for quick probability estimates and identifying outliers.

## 5. Key Properties

| Property | Formula | Meaning |
|----------|---------|---------|
| **Mean** | $E[X] = \mu$ | Center of the distribution |
| **Variance** | $\text{Var}(X) = \sigma^2$ | Spread around the mean |
| **Skewness** | 0 | Perfectly symmetric |
| **Shape** | Bell curve | Same shape regardless of μ and σ |

## 6. Central Limit Theorem

**Why the normal distribution is everywhere:**

The Central Limit Theorem states that the sum (or average) of many independent random variables approaches a normal distribution, regardless of the original distributions.

**Population genetics relevance:**
- Allele frequency changes over many generations approximate normality
- Quantitative traits influenced by many genes (polygenic) are normally distributed
- Sample means of genetic statistics become normally distributed with large sample sizes

**Example:** Even if individual mutations have non-normal effects, the total effect of many mutations on a trait tends toward normal.

## 7. Practical Example: Quantitative Trait

Suppose adult height in a population is normally distributed with μ = 170 cm and σ = 10 cm. What proportion of individuals are taller than 185 cm?

**Step 1: Standardize**

$$Z = \frac{185 - 170}{10} = 1.5$$

**Step 2: Look up or calculate**

Using a standard normal table or software, P(*Z* > 1.5) ≈ 0.067

**Result:** About 6.7% of the population is taller than 185 cm.

## 8. Connection to Other Distributions

The normal distribution emerges as a limit of several discrete distributions:

- **Binomial** → Normal when *n* is large (typically *n* > 30 and *np* > 5)
  - Approximation: *X* ~ N(*np*, *np*(1-*p*))
- **Poisson** → Normal when λ is large (typically λ > 20)
  - Approximation: *X* ~ N(λ, λ)

These approximations simplify calculations for large samples.

## 9. Intuitive Summary

The normal distribution is the "default" distribution for:
- Traits influenced by many small, independent factors
- Averages and sums of random variables
- Measurement error and biological variation

**Population genetics applications:**
- Modeling polygenic traits and quantitative genetics
- Approximating allele frequency distributions under drift
- Statistical inference (confidence intervals, hypothesis tests)
- Breeding values and heritability estimates

---

**Key insight:** While individual genetic events (mutations, Mendelian segregation) are discrete, their cumulative effects on populations and phenotypes often produce normally distributed outcomes.
