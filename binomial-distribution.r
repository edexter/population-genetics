################################################################################
# EXPLORING PROBABILITY DISTRIBUTIONS IN R
# A Beginner's Guide to Binomial, Poisson, and Normal Distributions
################################################################################

# HOW TO USE THIS SCRIPT:
# 1. Read through each section carefully
# 2. Run the code as-is first to see what happens
# 3. Then change the parameters (marked with "YOU CAN CHANGE")
# 4. Run it again and see how the plots change
# 5. Try to answer the questions at the end of each section

################################################################################
# PART 1: THE BINOMIAL DISTRIBUTION (COIN FLIPS)
################################################################################
# The binomial distribution models: "If I flip a coin N times, how many heads?"

# --- PARAMETERS YOU CAN CHANGE ---
number_of_flips <- 10         # How many coin flips in each experiment?
probability_heads <- 0.5      # What's the chance of heads? (0.5 = fair coin)
number_of_experiments <- 100  # How many experiments do we run?

# --- SIMULATE THE EXPERIMENTS ---
# rbinom() creates random binomial data
# We run "number_of_experiments" experiments
# Each experiment has "number_of_flips" flips
# Each flip has "probability_heads" chance of success (heads)
results <- rbinom(n = number_of_experiments, 
                  size = number_of_flips, 
                  prob = probability_heads)

# --- LOOK AT THE RAW RESULTS ---
print(results)

# --- MAKE A HISTOGRAM ---
# This shows us how many heads we got across all experiments
hist(results, 
     breaks = seq(-0.5, number_of_flips + 0.5, by = 1),  # One bar per possible outcome
     col = "lightblue",              # Color of bars
     main = paste("Binomial: Flipping", number_of_flips, "coins", number_of_experiments, "times"),
     xlab = "Number of Heads",       # X-axis label
     ylab = "Frequency (count)",     # Y-axis label
     probability = FALSE)            # Show counts, not probabilities

# --- ADD THE THEORETICAL CURVE ---
# This red line shows what the formula predicts
x_values <- 0:number_of_flips       # All possible outcomes (0 heads to max heads)

# dbinom() calculates the probability for each outcome using the formula (assuming a fair coin)
probabilities <- dbinom(x_values, size = number_of_flips, prob = 0.5)

# Scale it to match our histogram (convert probability to count)
scaled_probs <- probabilities * number_of_experiments

# Draw it as red points connected by lines
lines(x_values, scaled_probs, col = "red", lwd = 3, type = "b", pch = 19)

# --- PRINT SUMMARY STATISTICS ---
# Mean across all trials
mean(results)

# --- EXPLORING PROPORTIONS INSTEAD OF COUNTS ---
# Instead of counting heads, we can look at the PROPORTION of heads. This makes it
# it easier to compare experiment with different numbers of coin flips
proportions <- results / number_of_flips

# --- PRINT SUMMARY STATISTICS ---
# Mean across all trials
mean(proportions)

# Standard deviation across all trials
sd(proportions)


# Make a histogram of proportions
hist(proportions,
     col = "lightcoral",
     main = paste("Proportion of Heads (n =", number_of_flips, "flips)"),
     xlab = "Proportion of Heads",
     ylab = "Frequency (count)")

# QUESTIONS TO EXPLORE:
# 1. What happens when you change the number of trials to a very large or small number?
# 2. Model what happens when using an unfair coin.
# 3. What happens to the standard deviation if you increase or decrease the number of flips per trials?
# 4. Why do we look at the standard deviation of the proportion data instead of the raw number of heads per trial?

# --- REAL GENETICS EXAMPLE: PREDICTING OFFSPRING GENOTYPES ---
# When two heterozygotes (Aa × Aa) mate, what proportion of offspring are 
# homozygous recessive (aa)?
# From Mendelian genetics, we expect 1/4 (25%) to be aa

# Let's simulate this!
# Each offspring has a 0.25 probability of being aa (homozygous recessive)
probability_homozygous_recessive <- 0.25

# How many offspring should we count to get a good estimate?
number_of_offspring <- 20  # TRY CHANGING THIS: 10, 50, 100, 500

# Simulate ONE family with this many offspring
# We're asking: "Out of number_of_offspring kids, how many are aa?"
observed_aa <- rbinom(n = 1, 
                      size = number_of_offspring, 
                      prob = probability_homozygous_recessive)

# Convert to proportion
observed_proportion <- observed_aa / number_of_offspring

print("Expected proportion of aa offspring: 0.25 (25%)")
print(paste("Observed proportion:", round(observed_proportion, 3)))

# --- NOW SIMULATE MANY FAMILIES TO SEE THE PATTERN ---
# Let's imagine 100 different families, each with number_of_offspring kids
number_of_families <- 100

families_aa_counts <- rbinom(n = number_of_families,
                             size = number_of_offspring,
                             prob = probability_homozygous_recessive)

# Convert each family's count to a proportion
families_proportions <- families_aa_counts / number_of_offspring

# Plot the distribution of proportions across families
hist(families_proportions,
     col = "lavender",
     main = paste("Distribution of aa Proportions\n(", number_of_families, 
                  "families, each with", number_of_offspring, "offspring)"),
     xlab = "Proportion of aa Offspring",
     ylab = "Number of Families")

# Add a vertical line at the expected 0.25
abline(v = 0.25, col = "red", lwd = 3, lty = 2)

# Add a vertical line at the mean of our simulations
abline(v = mean(families_proportions), col = "blue", lwd = 2)

print(paste("Mean proportion across all families:", round(mean(families_proportions), 3)))

# GENETICS QUESTIONS TO EXPLORE:
# 1. Run this several times with number_of_offspring = 10
#    Notice how much the observed proportion jumps around!
# 2. Now change number_of_offspring to 100 and run it several times
#    Does it stay closer to 0.25?
# 3. Try number_of_offspring = 500. How tight is the distribution now?
# 4. KEY INSIGHT: Larger sample sizes give more reliable estimates!
#    This is why geneticists need large numbers of offspring to test ratios.
# 5. If you observed 10 aa out of 20 offspring (50%), is that really 
#    evidence against Mendelian inheritance? Or just sampling variation?
