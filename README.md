# FASTA-Toolkit: Genomic Sequence Analysis

A command-line Python utility designed to parse multi-FASTA files, detect Open Reading Frames (ORFs) across various reading frames, and track overlapping repeated nucleotide motifs.

*Note: This project was originally developed as part of an examination for the **Genomic Data Science** course by **Johns Hopkins University (JHU)** on Coursera. It has since been refactored, optimized, and structurally enhanced for inclusion in a professional engineering portfolio.*

---

## Core Features

* **FASTA File Parsing:** Efficiently streams and reconstructs multi-line DNA sequences from raw FASTA inputs into native key-value dictionary formats, safely managing trailing whitespaces or formatting quirks.
* **Sequence Metrics:** Computes detailed structural diagnostics, including exact total record counts, specific sequence lengths, and comprehensive tie-handling for the shortest and longest strings in a file.
* **Dynamic ORF Detection:** Systematically identifies complete Open Reading Frames (ORFs) across specific forward reading frames (1, 2, or 3) by enforcing standard molecular biology logic (strict `ATG` start to `TAA`/`TAG`/`TGA` stop translation boundaries).
* **Overlapping Motif Profiling:** Implements an $n$-length sliding window mechanism to track recurring sequence repeats on the forward strand, reverse-complement strand, and total combined orientations.

---

## What You Need to Do (Getting Started)

### 1. Prerequisites & Dependencies
Ensure you have Python 3.8+ installed on your system. This tool relies entirely on Python standard library modules, meaning zero external dependencies are strictly mandatory. However, a `requirements.txt` environment configuration is supplied to standardize portfolio deployments.

### 2. Installation
Clone the repository using your terminal or open the folder directly in VS Code:
### 3. Execution
Run the analysis script from your command line by passing your target FASTA dataset as an argument:
```
python fasta_analysis.py dna.example.fasta
```
### 4. Interactive Prompts
Upon launch, the script will interactively guide you through runtime customization:Select Reading Frame: Enter 1, 2, or 3 to establish the initial codon alignment step.Target ID Query: Input a specific FASTA sequence identifier (e.g., gi|142022655) to evaluate its isolated internal genomic features.Repeat Window Definition: Specify an integer sequence length $n$ (e.g., 12) to map structural repeats.