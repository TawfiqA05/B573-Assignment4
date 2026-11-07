# Assignment 4: Basic R Commands
# INFO-B 473/B573
# Tawfiq Abulail
#
# This script reads the chr1_GL383518v1_alt sequence I downloaded in Assignment 1.
# Part 1 prints letters from the sequence, Part 2 builds the reverse complement,
# Part 3 counts each letter per kilobase, and Part 4 puts the A/C/G/T counts in
# data frames and checks the row sums.
#
# How to run (from the folder that has the sequence file):
#   Rscript Assignment4.R
# or give the path to the file:
#   Rscript Assignment4.R path/to/chr1_GL383518v1_alt.fa
# In RStudio, set the working directory to this folder and click Source.
#
# Only base R is used, no packages needed.


########## Settings ##########

# file name from Assignment 1 (.fa or .fa.gz both work)
default_file <- "chr1_GL383518v1_alt.fa"
kb_size <- 1000
nucleotides <- c("A", "C", "G", "T")


########## Functions ##########

# Find the sequence file. Uses the command line path if one was given,
# otherwise looks for the .fa file and then the .fa.gz file.
find_sequence_file <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    candidates <- args[1]
  } else {
    candidates <- c(default_file, paste0(default_file, ".gz"))
  }
  for (f in candidates) {
    if (file.exists(f)) {
      return(f)
    }
  }
  stop("Could not find the sequence file. Looked for: ",
       paste(candidates, collapse = ", "),
       ". Put the file in this folder or pass its path to the script.")
}

# Read a FASTA file and return the whole sequence as one string.
# readLines() can also read a .gz file directly.
read_sequence <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) {
    stop("The file ", path, " is empty.")
  }

  headers <- grep("^>", lines)
  if (length(headers) == 0) {
    stop("No FASTA header line (starting with >) found in ", path)
  }
  if (length(headers) > 1) {
    cat("Note: file has", length(headers), "sequences, only the first one is used.\n")
    lines <- lines[1:(headers[2] - 1)]
  }

  # drop the header line and join the rest into one string
  seq_lines <- lines[-headers[1]]
  sequence <- paste(seq_lines, collapse = "")
  sequence <- gsub("[[:space:]]", "", sequence)

  if (nchar(sequence) == 0) {
    stop("The file has a header but no sequence.")
  }

  # check for letters that are not A, C, G, T or N (upper or lower case)
  bad <- unique(strsplit(gsub("[ACGTNacgtn]", "", sequence), "")[[1]])
  if (length(bad) > 0) {
    cat("Note: sequence has other letters:", paste(bad, collapse = " "), "\n")
  }
  return(sequence)
}

# Return the letter at one position (positions start at 1, like the assignment)
get_letter <- function(sequence, position) {
  get_range(sequence, position, position)
}

# Return the letters from start to end, including both ends
get_range <- function(sequence, start, end) {
  n <- nchar(sequence)
  if (!is.numeric(start) || !is.numeric(end) || start %% 1 != 0 || end %% 1 != 0) {
    stop("Positions must be whole numbers.")
  }
  if (start < 1 || end > n || start > end) {
    stop("Positions ", start, " to ", end, " are outside the sequence (length ", n, ").")
  }
  substr(sequence, start, end)
}

# Reverse the sequence and swap each base for its Watson-Crick-Franklin pair.
# Lower case letters stay lower case and N stays N.
reverse_complement <- function(sequence) {
  bases <- strsplit(sequence, "")[[1]]
  bases <- rev(bases)
  complement <- chartr("ACGTacgt", "TGCAtgca", paste(bases, collapse = ""))
  return(complement)
}

# Count how many times each letter shows up in a piece of sequence.
# The counts are kept in a named list, which works like a dictionary in R.
count_letters <- function(chunk) {
  counts <- list()
  for (letter in strsplit(chunk, "")[[1]]) {
    if (is.null(counts[[letter]])) {
      counts[[letter]] <- 1
    } else {
      counts[[letter]] <- counts[[letter]] + 1
    }
  }
  return(counts)
}

# Print helper: show a letter count dictionary in alphabetical order
sort_counts <- function(letter_counts) {
  counts <- unlist(letter_counts)
  counts[order(toupper(names(counts)), names(counts))]
}

# Turn one kilobase dictionary into a one row data frame with columns A, C, G, T.
# The file uses lower case for repeat regions, so "a" is added to "A" and so on.
kb_to_row <- function(letter_counts) {
  row <- c(A = 0, C = 0, G = 0, T = 0)
  for (base in nucleotides) {
    for (letter in c(base, tolower(base))) {
      if (!is.null(letter_counts[[letter]])) {
        row[base] <- row[base] + letter_counts[[letter]]
      }
    }
  }
  data.frame(A = row["A"], C = row["C"], G = row["G"], T = row["T"], row.names = NULL)
}


########## Main program ##########

main <- function() {

  ##### Part 1: read the sequence #####
  cat("Part 1\n")
  seq_file <- find_sequence_file()
  sequence <- read_sequence(seq_file)
  seq_length <- nchar(sequence)
  cat("File:", seq_file, "\n")
  cat("Sequence length:", seq_length, "bases\n")
  # the file keeps lower case letters (repeat regions), so some answers print in lower case

  # 1a and 1b
  cat("1a. 10th letter:", get_letter(sequence, 10), "\n")
  cat("1b. 758th letter:", get_letter(sequence, 758), "\n\n")


  ##### Part 2: reverse complement #####
  cat("Part 2\n")
  rev_comp <- reverse_complement(sequence)
  cat("Reverse complement length:", nchar(rev_comp), "bases\n")

  # 2a and 2b
  cat("2a. 79th letter:", get_letter(rev_comp, 79), "\n")
  cat("2b. 500th through 800th letters:\n")
  cat(get_range(rev_comp, 500, 800), "\n\n")


  ##### Part 3: letter counts for each kilobase #####
  # This uses the sequence that was already read in Part 1.
  cat("Part 3\n")
  n_kb <- ceiling(seq_length / kb_size)
  kb_counts <- list()

  for (i in 1:n_kb) {
    start <- (i - 1) * kb_size + 1
    end <- min(i * kb_size, seq_length)
    chunk <- substr(sequence, start, end)
    kb_counts[[paste0("kb_", i)]] <- count_letters(chunk)
  }

  cat("3a. Number of kilobases in the list:", length(kb_counts), "\n")
  cat("Letter counts for kb_1:\n")
  print(sort_counts(kb_counts[["kb_1"]]))
  cat("Letter counts for the last one, ", names(kb_counts)[n_kb], ":\n", sep = "")
  print(sort_counts(kb_counts[[n_kb]]))
  cat("Total letters counted:", sum(unlist(kb_counts)), "\n\n")


  ##### Part 4: data frames from the Part 3 list #####
  cat("Part 4\n")

  # 4a. data frame for the first 1000 base pairs
  first_kb <- kb_to_row(kb_counts[["kb_1"]])
  cat("4a. A, C, G, T counts in the first 1000 base pairs:\n")
  print(first_kb)

  # 4b. same thing for every kilobase in the dictionary
  kb_frames <- list()
  for (kb_name in names(kb_counts)) {
    kb_frames[[kb_name]] <- kb_to_row(kb_counts[[kb_name]])
  }
  cat("\n4b. Made", length(kb_frames), "data frames, one per kilobase. First three:\n")
  for (kb_name in names(kb_frames)[1:min(3, length(kb_frames))]) {
    cat(kb_name, "\n")
    print(kb_frames[[kb_name]])
  }

  # 4c. put each data frame from 4b in as its own row
  all_kb <- do.call(rbind, kb_frames)
  rownames(all_kb) <- names(kb_frames)
  cat("\n4c. Combined table:", nrow(all_kb), "rows and", ncol(all_kb), "columns\n")
  print(head(all_kb, 5))
  cat("...\n")
  print(tail(all_kb, 5))

  # 4d. sum of each row
  all_kb$Sum <- rowSums(all_kb[, nucleotides])
  cat("\n4d. Row sums (first and last five):\n")
  print(head(all_kb, 5))
  cat("...\n")
  print(tail(all_kb, 5))

  # save the full table so every row can be checked
  write.csv(all_kb, "kilobase_counts.csv")
  cat("Full table saved to kilobase_counts.csv\n")

  # 4e. compare the sums to the expected value
  expected <- kb_size
  off_rows <- all_kb[all_kb$Sum != expected, ]
  cat("\n4e. Expected sum for each row:", expected, "\n")
  if (nrow(off_rows) == 0) {
    cat("Every row adds up to", expected, "\n")
  } else {
    cat("Rows that do not equal ", expected, ":\n", sep = "")
    print(off_rows)
  }

  # 4e answers:
  # 1. Each list is one kilobase, so the expected sum is 1000.
  #
  # 2. Yes, one list. kb_1 through kb_182 all add up to 1000, but kb_183 only
  #    adds up to 439.
  #
  # 3. The sequence is 182,439 bases long, which is not a multiple of 1000.
  #    After 182 full kilobases there are only 439 bases left, so the last list
  #    is short. In general a sum can also come out lower than expected if the
  #    sequence has letters other than A, C, G and T, like N for unknown bases,
  #    because those are not counted in the 4 columns. This file has no N's.
  #    The file also writes repeat regions in lower case (a, c, g, t). If only
  #    upper case letters were counted, almost every row would be under 1000,
  #    so I added the lower case counts to the upper case ones.
}

# Run everything and stop with a clear message if something goes wrong
tryCatch(
  main(),
  error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    if (!interactive()) quit(status = 1)
  }
)
