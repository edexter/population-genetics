# Practical Exercise: Genetic Simulations with quantiNEMO

## Learning Objectives

By the end of this practical, you will be able to:
1. Install and run quantiNEMO on Linux systems
2. Create and modify simulation parameter files
3. Simulate mutation-drift equilibrium in populations
4. Model the effects of migration on genetic diversity
5. Interpret simulation outputs and compare to theoretical predictions
6. (Optional) Explore the effects of selection on allele frequencies

---

## Background

Population genetic theory makes predictions about how genetic variation changes over time due to mutation, drift, migration, and selection. However, these processes interact in complex ways that can be difficult to visualize or predict analytically.

**quantiNEMO** is a forward-time, individual-based simulation program that allows us to explicitly model these evolutionary processes. By running simulations, we can:
- Test theoretical predictions
- Explore parameter space (what happens when Ne is small vs. large?)
- Understand how different forces interact
- Generate expectations for comparison with real data

In this practical, we'll use quantiNEMO to explore fundamental population genetics concepts through simulation experiments.

---

## Section 1: Installation and Setup

### 1.1 System Requirements

**Operating System:**
- Linux (native or WSL2 for Windows users)
- macOS (with modifications to commands)

**Prerequisites:**

You'll need the `unzip` utility to extract the downloaded file.

**For Linux/WSL2:**
```bash
# Check if unzip is installed
which unzip

# If not installed, install it
sudo apt-get install unzip
```

**For macOS:**
```bash
# unzip is usually pre-installed on macOS
# If needed, install via Homebrew
brew install unzip
```

**Note:** quantiNEMO is pre-compiled and built with C++. No additional packages or compilation are required.

---


### 1.2 Download quantiNEMO

We'll follow the official installation instructions from the quantiNEMO website.

**Step 1: Navigate to your project directory**

```bash
# Navigate to your home directory
cd ~

# Create a directory for simulations
mkdir genetic_simulation
cd genetic_simulation

# Create subdirectories for organization
mkdir bin settings output results
```

**Step 2: Download quantiNEMO**

Visit the official download page: https://www2.unil.ch/popgen/softwares/quantinemo/getting_started.html

**For Linux users (including WSL2):**

```bash
# Navigate to the bin directory
cd bin

# Download the Linux version (it's a zip file)
wget https://www2.unil.ch/popgen/softwares/quantinemo/files/quantinemo_linux.zip

# Unzip the file (extracts to a directory called quantinemo_linux/)
unzip quantinemo_linux.zip

# This creates a directory with three files:
# - quantinemo (the executable)
# - quantinemo.ini (example settings file)
# - quantinemo.pdf (user manual)

# Move the executable to the bin directory
mv quantinemo_linux/quantinemo .

# Clean up the extracted directory and zip file
rm -r quantinemo_linux/
rm quantinemo_linux.zip

# Modify permissions so it can be launched
chmod +x quantinemo

# Verify the installation
ls -lh quantinemo
./quantinemo --version
```

**Expected output:**
- You should see the `quantinemo` executable (~5.5 MB)
- Running `./quantinemo --version` should display version information

---

### 1.3 Verify Installation (Health Check)

Let's test that quantiNEMO is working correctly.

**Test 1: Check version and help**

```bash
# Move back one level to the main simulation directory
cd ..

# Run quantiNEMO to see version information
./bin/quantinemo --version

# View help information
./bin/quantinemo --help
```

**Expected output:**
- Version information (e.g., "quantiNEMO version 2.x.x")
- List of available command-line options
- Information about settings files

**Test 2: Run a minimal test simulation**

Create a simple test settings file:

```bash
# Create a minimal test settings file
cat > settings/test.ini << 'EOF'
# Minimal quantiNEMO test simulation

#Number of generations
generations 1000

# metapopulation
patch_number 1
patch_capacity 1000  

# quantitative trait 
quanti_loci 1

# statistics 
stat {q.adlt.ho}
EOF

# View the file
cat settings/test.ini

# You can also interactively edit the file
nano settings/test.ini

```

**Run the test simulation:**

```bash
# Run quantiNEMO with the test settings
./bin/quantinemo settings/test.ini

# Check that it created output files
ls -lh
```

**Expected behavior:**
- quantiNEMO should run without errors
- You should see messages about simulation progress
- Output files may be created (depending on output settings)

**If you see errors:**
- Check that the executable has proper permissions (`ls -l bin/quantinemo`)
- Verify you're in the `genetic_simulation` directory
- Check that the settings file was created correctly (`cat settings/test.ini`)

---

### 1.4 Understanding quantiNEMO File Structure

quantiNEMO uses **settings files** (`.ini` format) to configure simulations. Let's understand the basic structure:

**Key components of a settings file:**

```ini
# Comments start with #

# Simulation parameters
generations 100          # How many generations to simulate
replicates 10           # How many independent simulations to run

# Population structure
patch_number 2          # Number of populations (demes)
patch_capacity 100      # Individuals per population

# Genetic architecture
quanti_loci 10          # Number of loci
quanti_alleles 2        # Alleles per locus (usually 2)
quanti_init_model 1     # How to initialize allele frequencies

# Evolutionary forces
mutation_rate 0.001     # Per-locus mutation rate
dispersal_rate 0.1      # Migration rate between populations

# Mating system
mating_system 1         # 1=random mating, 2=selfing, etc.

# Output options
stat_log_time {1}       # Output statistics at generation 1
stat {quanti.freq}      # Which statistics to output
```

**File organization in our project:**
```
genetic_simulation/
├── bin/                    # quantiNEMO executable
├── settings/               # .ini parameter files
│   ├── test.ini
│   ├── drift.ini
│   ├── migration.ini
│   └── selection.ini
├── output/                 # Raw simulation output
└── results/                # Processed results, plots
```

---

### 1.5 quantiNEMO Workflow Overview

**Typical workflow:**

1. **Design experiment** - What process do you want to explore?
2. **Create settings file** - Define parameters
3. **Run simulation** - Execute quantiNEMO
4. **Analyze output** - Process results (often with R)
5. **Visualize** - Create plots to interpret results
6. **Iterate** - Modify parameters and repeat

**Command structure:**

```bash
# Basic usage
./bin/quantinemo settings/my_settings.ini

# Redirect output to a file
./bin/quantinemo settings/my_settings.ini > output/my_output.txt

# Run in background (for long simulations)
nohup ./bin/quantinemo settings/my_settings.ini > output/my_output.txt 2>&1 &
```

---

### 1.6 Checkpoint: Are You Ready to Continue?

Before moving on to the simulation experiments, verify:

✅ quantiNEMO is installed in `~/genetic_simulation/bin/`
✅ You can run `./bin/quantinemo --version` successfully
✅ You ran the test simulation without errors
✅ You understand the basic structure of settings files
✅ Your directory structure is organized (bin/, settings/, output/, results/)

**If something isn't working:**
- Re-read the installation steps carefully
- Check file permissions (`ls -l bin/quantinemo`)
- Verify you're using the Linux version (not Windows .exe)
- Ask for help from your instructor

---

## Section 2: Experiment 1 - Mutation-Drift Equilibrium

# You will need to set up your ini file with the following lines:

cat > settings/drift_only.ini << 'EOF'
# General simulation settings
generations             10000       # Number of generations CHANGE THIS!
patch_capacity          100         # Population size CHANGE THIS!

# Neutral marker configuration
ntrl_loci               10          # Number of neutral loci
ntrl_all                10          # Number of alleles per locus
ntrl_ini_allele_model   0           # 0 = All loci start polymorphic, 1 = All loci started fixed,  
ntrl_mutation_rate      0.0         # Mutation rate

# Statistics output
stat             {n.adlt.nbAll_p1}  # Mean alleles/locus
stat_log_time           1           # Record statistics every generation
EOF

# Output file naming
filename                drift_only

In this experiment we are going to explore how the processes of drift in mutation act in opposite directions to either increase or decrease genetic variation in a population. For this experiment we will make some simplifying assumptions:

1. Constant population size. The total number of individuals does not change from one generation to the next. This is a common assumption in population genetic models.

2. Neutral markers only. The genetic markers that we simulate are not subject to natural selection, meaning different alleles all produce the same phenotype. Often when trying to understand population processes we seek to use neutral markers because selection can distort the expected signals of neutral processes

3. Indepedent markers. There is no linkage disequilbirum between markers and thus each one is independent

4. No migration. New genetic variants can only arise in the population from mutation. In future simulations we will relax this assumption and allow for migration.

Exercise 1. First we explore the effect of drift alone. Start with a single population of size n with 10 neutral markers, each with 10 alleles. Set the mutation rate to 0. Run the simulation for n generations and monitor how the average number of alleles changes over time. Experiment with different values to answer the following questions:

1. What effect does drift have on genetic diversity?
2. How does drift interact with population size?
3. If you could run your simulation for an infinitly large number of generations, what value would you expect the average number of alleles per locus to approach?
3. What is the minimum size this population needs to be to maintain an average of 2 alleles at each locus for more than 5000 generations? Demonstrate your answer with a single R plot.

Exercise 2. Now we introduce the force of mutation to counterbalance drift. Create a new ini file for quantiNemo named "drift_mutation.ini". Explore the effect of different mutation rates on your simulation and answer the following questions:

cat > settings/drift_mutation.ini << 'EOF'
# General simulation settings
generations             10000       # Number of generations CHANGE THIS!
patch_capacity          100         # Population size CHANGE THIS!

# Neutral marker configuration
ntrl_loci               10          # Number of neutral loci
ntrl_all                10          # Number of alleles per locus
ntrl_ini_allele_model   0           # 0 = All loci start polymorphic, 1 = All loci started fixed,  
ntrl_mutation_rate      0.0         # Mutation rate

# Statistics output
stat             {n.adlt.nbAll_p1}  # Mean alleles/locus
stat_log_time           1           # Record statistics every generation
EOF

# Output file naming
filename                drift_mutation
EOF

1. What effect does mutation have on genetic diversity?
2. For the population size that you found in answering questions 3 above, how long does it take for the population to reach mutation-drift equlibrium? What is the value of genetic diversity it settles at (alleles/locus)
3. Does the population approach mutation drift equilibrium at a constant speed, or does it vary depending on how far the current value is from the equilibrium value? Does your answer change depending on whether you start the simulation with all loci polymorphic or all loci monomorphic?

Exercise 3. Now we introduce a third force into our scenario: Migration. In this scenario we will start with 2 populations of the same size, mutation rate, and initial genetic configuration (10 loci with 10 alleles all maximal polymorphic). We will have a migration rate, but the migration rate will be set to 0 for our first simulation, so that the populations will evolve in parallel in a completely independent matter. At each generation in the simulation we will calculate the statistic Fst which measures how much population structure exists in the metapopulation and how diverged these two populations are.

cat > settings/drift_mutation_migration.ini << 'EOF'
# General simulation settings
generations             5000       # Number of generations CHANGE THIS!
patch_capacity          1000         # Population size CHANGE THIS!
patch_number            2

# Neutral marker configuration
ntrl_loci               100          # Number of neutral loci
ntrl_all                2          # Number of alleles per locus
ntrl_ini_allele_model   0           # 0 = All loci start polymorphic, 1 = All loci started fixed,  
ntrl_mutation_rate      0.001         # Mutation rate

# Migration configuration
dispersal_rate 0.001

# Statistics output
stat             {n.adlt.nbAll_p1 n.adlt.nbAll_p2 n.adlt.fst}  # Mean alleles/locus
stat_log_time           1           # Record statistics every generation

# Output file naming
filename                drift_mutation_migration
EOF

1. Plot how FST evolves through time in your simulation. Does it at some point reach a plateau? Is that biologically reasonable or an artifact of the simulation?

2. Does the change in FST (and hence differentiation) happen more quickly or more slowly with larger or smaller populations?

Exercise 4. We saw that isolation will cause populations to diverge across time, now we ask the question of what happens when we connect them via migration. Experiment with increasing values of the migration parameter to asnwer the following questions

1. How does migration affect the rate at which populations drift apart from each other?

2. How does migration change the level of diversity that the populations can maintain. Do the number of alleles maintained in a population scale with the migration rate?


---

## Additional Resources

### Documentation
- **quantiNEMO website:** https://www2.unil.ch/popgen/softwares/quantinemo/
- **User manual:** Available on the website (comprehensive guide to all parameters)
- **Tutorials:** R markdown tutorials provided by developers

### Key Papers
- Neuenschwander et al. (2019). QuantiNemo 2: a Swiss knife to simulate complex demographic and genetic scenarios. *Bioinformatics*, 35(5), 886-888.

---

**End of Section 1**

**Next:** We'll explore mutation-drift equilibrium by simulating populations of different sizes and comparing observed heterozygosity to theoretical predictions.
