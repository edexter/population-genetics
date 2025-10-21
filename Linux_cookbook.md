# Linux Cookbook for Bioinformatics

A quick reference guide for essential Linux commands you'll need in this course. Commands are organized by complexity and frequency of use:

- **Basic Commands** - You will use these constantly. **Memorize these ASAP because** stopping to look them up frequently will significantly slow down your bioinformatics work. You will need to write out basic linux commands by memory for the final examination.
- **Intermediate Commands** - You will use some of these regularly, but you will not be required to memorize them for the course. As long as you know these commands exist, when to use them, and where to look up how to construct them, you're fine.
- **Advanced Commands** - These involve complex syntax that you'll need to carefully construct with help from a reference document. Just be aware that they exist and are available to you.

---

## Table of Contents

1. [Introduction to Linux](#introduction-to-linux)
2. [Basic Commands](#basic-commands)
   - [Navigating directories (cd)](#navigating-directories-cd)
   - [Listing and creating directories (ls, mkdir)](#listing-and-creating-directories-ls-mkdir)
   - [Viewing file contents (less, cat, head, tail)](#viewing-file-contents-less-cat-head-tail)
   - [Counting lines and words (wc)](#counting-lines-and-words-wc)
   - [Editing files (nano)](#editing-files-nano)
   - [Moving/renaming files (mv)](#movingrenaming-files-mv)
   - [Copying files (cp)](#copying-files-cp)
   - [Deleting files (rm)](#deleting-files-rm)
   - [Using wildcards (*)](#using-wildcards-)
   - [Making scripts executable (chmod)](#making-scripts-executable-chmod)
   - [Redirecting output (>, >>)](#redirecting-output--)
   - [Piping commands (|)](#piping-commands-)
3. [Advanced Linux](#advanced-linux)
   - [Searching file contents (grep)](#searching-file-contents-grep)
   - [Stream editing (sed, awk)](#stream-editing-sed-awk)
4. [Intermediate Commands](#intermediate-commands)
   - [Installing software](#installing-software)
   - [Connecting to servers (ssh)](#connecting-to-servers-ssh)
   - [Running background jobs (screen)](#running-background-jobs-screen)

---

## Introduction to Linux

### What is Linux and why do we use it?

Linux is a Unix-based operating system that dominates scientific computing and bioinformatics. Most bioinformatics tools are designed for Linux, high-performance computing clusters run Linux, and the command-line interface provides powerful tools for processing large datasets efficiently.

### Installing and configuring Linux on your local machine

**For practice during this course:**
- Use [JSLinux](https://bellard.org/jslinux/) in your browser (no installation needed)

**For your own computer:**
- **Windows users:** Install [WSL2 (Windows Subsystem for Linux)](https://docs.microsoft.com/en-us/windows/wsl/install)
- **Mac users:** Use the built-in Terminal application (macOS is Unix-based)
- **Linux users:** You're already set up!

---

## Basic Commands

### Navigating directories (cd)

Change your current working directory.

**Syntax:**
```bash
cd [directory_path]
```

**Examples:**
```bash
cd /home/username/data          # Go to absolute path
cd projects                     # Go to subdirectory in current location
cd ..                           # Go up one directory level
cd ../..                        # Go up two directory levels
cd ~                            # Go to your home directory
cd -                            # Go back to previous directory
```

**Related command:**
```bash
pwd                             # Print working directory (show where you are)
```

---

### Listing and creating directories (ls, mkdir)

List directory contents and create new directories to organize your work.

**Syntax:**
```bash
ls [options] [directory]
mkdir [directory_name]
```

**Examples:**
```bash
ls                              # List contents of current directory
ls -lh                          # List with details (size, permissions, date)
ls -a                           # List all files including hidden ones
ls -lht                         # List sorted by modification time

mkdir analysis                  # Create new directory
mkdir -p project/data/raw       # Create nested directories
```

---

### Viewing file contents (less, cat, head, tail)

Display the contents of text files without editing them.

**Syntax:**
```bash
less [filename]                 # View file interactively (recommended for large files)
cat [filename]                  # Print entire file to screen
head [filename]                 # Show first 10 lines
tail [filename]                 # Show last 10 lines
```

**Examples:**
```bash
less sequences.fasta            # View FASTA file (press 'q' to quit, space to scroll)
cat results.txt                 # Display entire file
cat sample1.txt sample2.txt     # Display multiple files concatenated

head sequences.fasta            # Check first 10 lines (verify file format)
head -n 20 data.txt             # Show first 20 lines
tail -n 50 log.txt              # Show last 50 lines (check job status)
tail -f analysis.log            # Follow file as it grows (monitor running jobs)
```

**Navigation in `less`:**
- `Space` or `f` - forward one page
- `b` - backward one page
- `/search_term` - search forward
- `q` - quit

---

### Counting lines and words (wc)

Count lines, words, and characters in files. Especially useful for counting sequences or checking file sizes.

**Syntax:**
```bash
wc [options] [filename]
```

**Examples:**
```bash
wc -l sequences.fasta           # Count lines in file
wc -l *.fastq                   # Count lines in all FASTQ files
wc -c large_file.txt            # Count characters (file size)
```

**Common options:**
- `-l` - Count lines only
- `-w` - Count words only
- `-c` - Count characters (bytes)

---

### Editing files (nano)

Simple text editor for creating and modifying files from the command line.

**Syntax:**
```bash
nano [filename]
```

**Examples:**
```bash
nano analysis_notes.txt         # Create or edit a file
nano config.txt                 # Edit configuration file
```

**Key commands in nano:**
- `Ctrl+O` - Save (write out)
- `Ctrl+X` - Exit
- `Ctrl+K` - Cut line
- `Ctrl+U` - Paste
- `Ctrl+W` - Search

---

### Moving/renaming files (mv)

Move files to different locations or rename them.

**Syntax:**
```bash
mv [source] [destination]
```

**Examples:**
```bash
mv old_name.txt new_name.txt    # Rename a file
mv data.txt results/            # Move file to directory
mv sample*.fasta raw_data/      # Move multiple files using wildcards
mv file.txt ../                 # Move file up one directory
```

---

### Copying files (cp)

**Syntax:**
```bash
cp [source] [destination]
```

**Examples:**
```bash
cp original.txt copy.txt        # Copy single file
cp data.txt backup/             # Copy file to directory
cp -r my_folder my_folder_backup  # Copy directory and all contents
cp *.fastq analysis/            # Copy multiple files using wildcards
```

**Common options:**
- `-r` - Recursive (required for copying directories)
- `-i` - Interactive (ask before overwriting)

---

### Deleting files (rm)

Remove files and directories permanently (cannot be undone).

**Syntax:**
```bash
rm [filename]
```

**Examples:**
```bash
rm unwanted.txt                 # Delete single file
rm file1.txt file2.txt          # Delete multiple files
rm *.tmp                        # Delete all files matching pattern
rm -r old_analysis/             # Delete directory and contents
rm -i important_*               # Delete with confirmation prompt
```

**⚠️ WARNING:**
- There is no "recycle bin" in Linux - deleted files are gone forever
- Be extremely careful with `rm -r` and wildcards
- Always double-check before pressing Enter

**Common options:**
- `-r` - Recursive (required for directories)
- `-i` - Interactive (ask confirmation for each file)
- `-f` - Force (no confirmation, use with extreme caution)

---

### Using wildcards (*)

Wildcards let you match multiple files at once, making operations much faster but also riskier.

**Common wildcards:**
- `*` - Matches any number of characters
- `?` - Matches exactly one character
- `[abc]` - Matches any one of the enclosed characters

**Examples:**
```bash
ls *.fasta                      # List all FASTA files
ls sample_?.txt                 # List sample_1.txt, sample_A.txt, etc.
ls sample_[1-5].txt             # List sample_1.txt through sample_5.txt

cp *.fastq raw_reads/           # Copy all FASTQ files
rm temp_*                       # Delete all files starting with "temp_"

cat sample*.txt > combined.txt  # Concatenate all matching files
```

**⚠️ WARNING:**
- Always test wildcards with `ls` before using them with `rm` or `mv`
- Example: `ls *.tmp` first, then `rm *.tmp` if the list looks correct

---

### Making scripts executable (chmod)

Change file permissions to make scripts runnable. You'll need this when writing your own analysis scripts.

**Syntax:**
```bash
chmod [permissions] [filename]
```

**Examples:**
```bash
chmod +x my_script.sh           # Make script executable
chmod +x analysis.py            # Make Python script executable
./my_script.sh                  # Run the script

chmod 755 script.sh             # Make executable, readable by all
chmod 644 data.txt              # Make readable/writable by owner, readable by others
```

**Common permissions:**
- `+x` - Add execute permission
- `755` - Owner can read/write/execute, others can read/execute
- `644` - Owner can read/write, others can only read

---

### Redirecting output (>, >>)

Save command output to files instead of displaying on screen. Essential for saving analysis results.

**Syntax:**
```bash
command > output.txt            # Write output to file (overwrite)
command >> output.txt           # Append output to file
```

**Examples:**
```bash
ls -l > file_list.txt           # Save directory listing to file
echo "Analysis complete" > status.txt

grep ">" sequences.fasta > headers.txt        # Save FASTA headers
cat file1.txt file2.txt > combined.txt        # Combine files

echo "Sample 1 done" >> log.txt               # Append to log file
date >> log.txt                               # Add timestamp to log
```

**⚠️ WARNING:**
- `>` overwrites the file completely
- `>>` adds to the end of the file
- Double-check which one you need!

---

### Piping commands (|)

Send the output of one command directly to another command. This is the foundation of powerful one-line analyses in Linux.

**Syntax:**
```bash
command1 | command2
```

**Examples:**
```bash
ls -l | less                    # View long directory listing page by page
cat large_file.txt | head       # Show first 10 lines of file

ls | wc -l                      # Count number of files in directory
cat sequences.fasta | grep ">" | wc -l    # Count sequences in FASTA file

history | grep "grep"           # Search your command history
```

**Why this matters:**
Piping lets you build complex analyses from simple building blocks without creating intermediate files.

---

## Advanced Linux

This section covers powerful tools for text processing and data manipulation. These commands can be complex, but knowing they exist will help you find solutions when you need them.

### Searching file contents (grep)

Search for patterns in files. Extremely useful for filtering data, finding sequences, or extracting specific lines.

**Syntax:**
```bash
grep [pattern] [filename]
```

**Examples:**
```bash
grep "ATCG" sequences.fasta     # Find lines containing "ATCG"
grep ">" sequences.fasta        # Find all FASTA headers
grep -c ">" sequences.fasta     # Count how many headers (sequences)
grep -v ">" sequences.fasta     # Show all lines NOT containing ">" (sequences only)

grep "error" logfile.txt        # Find errors in log
grep -i "error" logfile.txt     # Case-insensitive search
grep -n "chr1" variants.vcf     # Show line numbers with matches
```

**Common options:**
- `-i` - Case-insensitive search
- `-v` - Invert match (show lines that DON'T match)
- `-c` - Count matching lines
- `-n` - Show line numbers

**Regular expressions (basic patterns):**
```bash
grep "^>" sequences.fasta       # Lines starting with ">"
grep "ATCG$" sequences.fasta    # Lines ending with "ATCG"
grep "A.G" sequences.fasta      # "A", any character, then "G"
grep "AT[CG]G" sequences.fasta  # "AT" followed by C or G, then "G"
```

`grep` is one of the most-used commands in bioinformatics. It's worth learning well.

---

### Stream editing (sed, awk)

`sed` and `awk` are powerful tools for transforming and processing text files. They can perform complex operations on data streams without loading entire files into memory - crucial for large bioinformatics files.

**When to use these:**
- Reformatting data files
- Extracting specific columns
- Find-and-replace operations
- Computing statistics on tabular data

**⚠️ Note:** These tools have a steep learning curve. The examples below show what's possible, but you'll likely need to look up syntax when you need them.

#### sed (Stream Editor)

Used for find-and-replace and simple text transformations.

**Examples:**
```bash
sed 's/old/new/' file.txt       # Replace first occurrence per line
sed 's/old/new/g' file.txt      # Replace all occurrences
sed -i 's/old/new/g' file.txt   # Edit file in-place

sed -n '1,10p' file.txt         # Print lines 1-10
sed '/pattern/d' file.txt       # Delete lines matching pattern
```

#### awk (Pattern Processing Language)

Used for column-based operations and more complex text processing.

**Examples:**
```bash
awk '{print $1}' data.txt       # Print first column
awk '{print $1, $3}' data.txt   # Print columns 1 and 3
awk '$3 > 100' data.txt         # Print lines where column 3 > 100

# Calculate average of column 2
awk '{sum+=$2} END {print sum/NR}' data.txt

# Filter VCF file for quality > 30
awk '$6 > 30' variants.vcf
```

**Why these matter:**
In bioinformatics, you'll frequently need to:
- Extract columns from tab-delimited files
- Filter based on numeric thresholds
- Reformat data between tools
- Process files too large for Excel

When you encounter these needs, remember that `awk` and `sed` exist, then search for specific examples online.

---

## Intermediate Commands

### Installing software

Most bioinformatics tools can be installed via package managers or conda. Package managers handle dependencies automatically.

**Using apt (Ubuntu/Debian):**
```bash
sudo apt update                 # Update package list
sudo apt install samtools       # Install samtools
```

**Using conda (recommended for bioinformatics):**

First install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), then:

```bash
conda install -c bioconda bwa          # Install BWA aligner
conda install -c bioconda fastqc       # Install FastQC
conda create -n myproject python=3.9   # Create isolated environment
conda activate myproject               # Activate environment
```

**From source (when necessary):**
```bash
wget https://example.com/tool.tar.gz   # Download
tar -xzf tool.tar.gz                   # Extract
cd tool/
./configure                            # Configure
make                                   # Compile
sudo make install                      # Install
```

---

### Connecting to servers (ssh)

SSH (Secure Shell) lets you connect to remote Linux servers securely, essential for accessing computing clusters.

**Syntax:**
```bash
ssh username@hostname
```

**Examples:**
```bash
ssh jsmith@cluster.university.edu       # Connect to remote server
ssh -p 2222 user@server.com            # Connect using non-standard port
```

**Copying files to/from remote servers:**
```bash
scp local_file.txt user@server:~/          # Copy TO server
scp user@server:~/results.txt ./           # Copy FROM server
scp -r my_folder/ user@server:~/projects/  # Copy directory recursively
```

**Tips:**
- You'll be prompted for password (or use SSH keys for passwordless login)
- Type `exit` or press `Ctrl+D` to disconnect
- Use `scp` (secure copy) to transfer files between machines

---

### Running background jobs (screen)

Screen allows you to run long-running analyses that continue even if you disconnect from the server.

**Basic usage:**
```bash
screen                          # Start new screen session
screen -S analysis              # Start named session
```

**Inside a screen session:**
```bash
# Run your long analysis here
bwa mem reference.fa reads.fastq > alignment.sam
```

**Key commands:**
- `Ctrl+A` then `D` - Detach from session (job keeps running)
- `screen -r` - Reattach to session
- `screen -ls` - List all sessions
- `screen -r analysis` - Reattach to specific named session
- `exit` - End the screen session (while inside it)

**Example workflow:**
```bash
# Start analysis
ssh user@cluster.edu
screen -S mapping
bwa mem ref.fa reads.fq > output.sam   # This will take hours
# Press Ctrl+A then D to detach
exit                                    # Disconnect from server

# Later, check on it
ssh user@cluster.edu
screen -r mapping                       # Reattach to see progress
```

---

## Additional Resources

- **Practice environment:** [JSLinux](https://bellard.org/jslinux/)
- **Command reference:** `man [command]` (e.g., `man ls` shows the manual for ls)
- **Get help:** Most commands support `--help` flag (e.g., `ls --help`)

---

## Quick Reference Card

| Task | Command |
|------|---------|
| Where am I? | `pwd` |
| List files | `ls` or `ls -lh` |
| Change directory | `cd [directory]` |
| Go home | `cd ~` |
| Go up one level | `cd ..` |
| View file | `less [file]` |
| Edit file | `nano [file]` |
| Copy file | `cp [source] [dest]` |
| Move/rename | `mv [source] [dest]` |
| Delete file | `rm [file]` |
| Delete directory | `rm -r [directory]` |
| Get help | `man [command]` or `[command] --help` |
