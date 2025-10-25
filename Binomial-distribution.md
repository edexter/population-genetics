The Binomial Distribution — A Primer for Population Genetics
1. Concept Overview
The binomial distribution models the probability of obtaining a specific number of successes in a fixed number of independent trials, where each trial has only two possible outcomes:

Success (e.g., getting heads in a coin toss, or inheriting a particular allele)
Failure (e.g., getting tails, or not inheriting that allele)

Key assumptions:

Each trial is independent (the outcome of one doesn't affect the others)
Each trial has the same probability of success, p
The number of trials, n, is fixed in advance

2. The Formula
The probability of getting exactly k successes in n trials is:
P(X=k)=(nk)pk(1−p)n−kP(X = k) = \binom{n}{k} p^k (1 - p)^{n - k}P(X=k)=(kn​)pk(1−p)n−k
Where:

X = the random variable representing the number of successes
n = total number of trials
k = the specific number of successes we're calculating the probability for
p = probability of success on each trial
(1−p)(1 - p)
(1−p) = probability of failure on each trial

(nk)=n!k!(n−k)!\binom{n}{k} = \frac{n!}{k!(n-k)!}
(kn​)=k!(n−k)!n!​ = the binomial coefficient, representing the number of ways to choose which *k* trials out of *n* are successes


3. Example: Flipping a Fair Coin 5 Times
Suppose you flip a fair coin (p = 0.5) 5 times. What is the probability of getting exactly 3 heads?
P(X=3)=(53)(0.5)3(0.5)2P(X = 3) = \binom{5}{3} (0.5)^3 (0.5)^{2}P(X=3)=(35​)(0.5)3(0.5)2
Step 1: Calculate the binomial coefficient
(53)=5!3!(5−3)!=5×4×3!3!×2!=5×42×1=10\binom{5}{3} = \frac{5!}{3!(5-3)!} = \frac{5 \times 4 \times 3!}{3! \times 2!} = \frac{5 \times 4}{2 \times 1} = 10(35​)=3!(5−3)!5!​=3!×2!5×4×3!​=2×15×4​=10
There are 10 different ways to get exactly 3 heads in 5 tosses.
Step 2: Plug into the formula
P(X=3)=10×(0.5)3×(0.5)2=10×(0.5)5=10×132=1032=0.3125P(X = 3) = 10 \times (0.5)^3 \times (0.5)^2 = 10 \times (0.5)^5 = 10 \times \frac{1}{32} = \frac{10}{32} = 0.3125P(X=3)=10×(0.5)3×(0.5)2=10×(0.5)5=10×321​=3210​=0.3125
Result: There's about a 31.25% chance of getting exactly 3 heads in 5 flips.
4. Understanding "5 Choose 3" — The Binomial Coefficient
The term (53)\binom{5}{3}
(35​), read as "5 choose 3," represents the number of different ways to select 3 items from 5, regardless of order. In our context, it's the number of different sequences of 5 coin flips that contain exactly 3 heads.

All 10 possible sequences with exactly 3 heads:

H H H T T
H H T H T
H H T T H
H T H H T
H T H T H
H T T H H
T H H H T
T H H T H
T H T H H
T T H H H

That's 10 total sequences — exactly what (53)=10\binom{5}{3} = 10
(35​)=10 tells us.

5. Key Properties
PropertyFormulaMeaningMeanE[X]=npE[X] = np
E[X]=npAverage number of successes across many repetitionsVarianceVar(X)=np(1−p)\text{Var}(X) = np(1 - p)
Var(X)=np(1−p)Spread of possible outcomesShapeSymmetric when p = 0.5Skewed right when p < 0.5; skewed left when p > 0.5
6. Intuitive Summary
The binomial distribution tells us how likely it is to get a certain number of successes when we repeat a binary (yes/no) experiment a fixed number of times.
Example insights for 5 coin flips:

Getting 0 or 5 heads is relatively rare (probability ≈ 3% each)
Getting 2 or 3 heads is most likely (probability ≈ 31% each)
All possible outcomes (0, 1, 2, 3, 4, or 5 heads) have probabilities that sum to 1

Connection to population genetics: The binomial distribution is fundamental for modeling genetic drift, allele sampling, and many other stochastic processes in population genetics where we track the number of copies of an allele across discrete generations or individuals.
