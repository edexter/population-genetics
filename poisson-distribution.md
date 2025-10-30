# The Poisson Distribution — A Primer for Population Genetics

## 1. Concept Overview

The Poisson distribution models the probability of observing a specific number of rare, independent events occurring in a fixed interval of time or space, when we know the average rate at which these events occur.

**Key characteristics:**
- Events occur independently at a constant average rate (λ)
- Two events cannot occur at exactly the same instant

**Population genetics examples:**
- Number of mutations occurring in a DNA sequence per generation
- Number of new alleles arising in a population over time
- Number of coalescent events in a small time interval

## 2. The Formula

The probability of observing exactly *k* events when the average rate is λ is:

$$P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}$$

**Where:**
* *X* = the random variable representing the number of events
* *k* = the specific number of events (0, 1, 2, 3, ...)
* λ (lambda) = the expected number of events (the mean)
* *e* ≈ 2.71828 (Euler's number)
* *k*! = k factorial

## 3. Example: Mutations in a DNA Sequence

Suppose mutations occur at an average rate of λ = 2.5 mutations per gene per million years. What is the probability of observing exactly 3 mutations?

$$P(X = 3) = \frac{(2.5)^3 e^{-2.5}}{3!}$$

**Calculation:**

$$P(X = 3) = \frac{15.625 \times 0.0821}{6} \approx 0.214$$

**Result:** About a 21.4% chance of observing exactly 3 mutations.

## 4. Connection to the Binomial Distribution

The Poisson distribution is a special case of the binomial when:
- *n* (number of trials) is very large
- *p* (probability per trial) is very small
- The product *np* = λ remains moderate

**Example:** A gene of 10,000 base pairs with mutation rate 0.00025 per base pair gives λ = 2.5 expected mutations. The Poisson is much simpler to calculate than the binomial here!

## 5. Key Properties

| Property | Formula | Meaning |
|----------|---------|---------|
| **Mean** | $E[X] = \lambda$ | Average number of events |
| **Variance** | $\text{Var}(X) = \lambda$ | Equals the mean! |
| **Shape** | Skewed right when λ is small | Approaches normal as λ increases |

**Important:** Mean equals variance is a distinctive property that can test whether data follow a Poisson process.

## 6. Practical Example: Testing Neutrality

Under neutral evolution, synonymous mutations should accumulate following a Poisson process. If we expect λ = 10 synonymous mutations but observe only *k* = 3:

$$P(X = 3) = \frac{10^3 e^{-10}}{3!} \approx 0.0076$$

This low probability (< 1%) might suggest selection or other non-neutral processes.

## 7. Intuitive Summary

The Poisson distribution answers: **"If rare events happen at a constant average rate, what's the probability of seeing exactly *k* events?"**

**Population genetics applications:**
- Modeling mutational processes and molecular evolution
- Describing polymorphism distributions under neutral theory
- Analyzing rare genetic events like duplications or rearrangements

---

**Quick comparison:** Use Poisson instead of binomial when *n* > 100, *p* < 0.01, and you're modeling rare events with known average rate λ.
