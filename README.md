# Assignment 4: Basic R Commands

Name: Tawfiq Abulail
Course: INFO-B 473/B573
Language: R (written and tested with R version 4.3.3)
Date submitted: November 7, 2026

## Description

This program reads the chr1_GL383518v1_alt sequence from Assignment 1 and works through the four parts of the assignment.

Part 1 reads the full sequence from the FASTA file and prints the 10th and 758th letters.

Part 2 makes the reverse complement of the sequence (reverse it, then swap A with T and C with G) and prints the 79th letter and the 500th through 800th letters.

Part 3 splits the sequence into kilobases and builds a list (a dictionary in R) that holds how many times each letter shows up in each kilobase.

Part 4 reads that list and makes a data frame of A, C, G and T counts for the first 1000 base pairs, does the same for every kilobase, stacks them into one table with one row per kilobase, and adds up each row. The answers to the 4e questions are in the comments at the end of the main program.

The file writes repeat regions in lower case, so lower case letters are kept in Parts 1 and 2 and counted together with the upper case letters in Part 4.

## Files

Assignment4.R: the R script for all four parts
chr1_GL383518v1_alt.fa: the sequence file downloaded in Assignment 1 (182,439 bases)
README.md: this file

## Libraries

None. Only base R is used.

## How to run

1. Put Assignment4.R and chr1_GL383518v1_alt.fa in the same folder.
2. From a terminal in that folder run: Rscript Assignment4.R
3. To use a file somewhere else, pass the path: Rscript Assignment4.R path/to/chr1_GL383518v1_alt.fa
4. In RStudio, set the working directory to the folder and click Source.

The script also works with the compressed chr1_GL383518v1_alt.fa.gz file. If the file is missing, empty or the positions are outside the sequence, it prints an error message and stops.

## Output

The answers for each part print to the screen.
kilobase_counts.csv: the full table from Part 4 with the A, C, G, T counts and the sum for all 183 kilobases.

## Results

1a. 10th letter: c
1b. 758th letter: a
2a. 79th letter of the reverse complement: G
2b. 500th through 800th letters of the reverse complement are printed by the script (301 letters).
The sequence has 183 kilobases. kb_1 through kb_182 each add up to 1000. kb_183 adds up to 439 because the sequence length is not a multiple of 1000.
