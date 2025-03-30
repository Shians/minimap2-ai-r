# minimap2 R Package

This is a NON-FUNCTIONAL EXPERIMENTAL PACKAGE created using AI to implement an R interface to minimap2.

## Installation

To install the package, you'll need to have R and Rcpp installed. You can install the package using:

```R
install.packages("devtools")
devtools::install_github("shians/minimap2-ai-r")
```

## Usage

Here's a basic example of how to use the package:

```R
library(minimap2)

# Create a minimap2 aligner object
aligner <- Minimap2Aligner(
    reference_file = "path/to/reference.fa",
    preset = "lr:hq",  # for short reads
    n_threads = 4
)

# Align a sequence
results <- align_sequences(
    aligner,
    query_seq = "ATCGATCGATCG",
    query_name = "read1"
)

# Clean up
destroy_aligner(aligner)
```

