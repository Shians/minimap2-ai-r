# minimap2 R Package

This is a NON-FUNCTIONAL EXPERIMENTAL PACKAGE created using AI to implement an R interface to minimap2.

## Installation

To install the package, you'll need to have R and Rcpp installed. You can install the package using:

```R
install.packages("remotes")
remotes::install_github("shians/minimap2-ai-r")
```

## Usage

A binary of minimap2 is now available at

```
system.file("bin/minimap2", package = "minimap2", mustWork = TRUE)
```

Here's a basic example of how to use the R interface of the package, provide your own human reference genome if you have one:

```R
library(minimap2)
library(tibble)

reads <- readr::read_lines(
    file = system.file("extdata/test.fastq.gz", package = "minimap2", mustWork = TRUE),
    n_max = 4000
)

reads_df <- tibble(
    name = reads[seq(1, length(reads), by = 4)],
    seq = reads[seq(2, length(reads), by = 4)],
    quals = reads[seq(4, length(reads), by = 4)]
)

genome_path <- tempfile()

# Might need more time to download genome
options(timeout = max(600, getOption("timeout")))
download.file(
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz",
    destfile = genome_path
)

aligner <- build_index(
    genome_path,
    preset = "map-ont",
    n_threads = 8
)

x <- aligner_map(
    aligner,
    query_seq = reads_df$seq,
    query_name = reads_df$name,
    query_quals = reads_df$quals
)

x
#> AlignmentRecords object
#> Header lines: 18 
#> Total records: 1720 
#> Mapped records: 261 
#> Unmapped records: 1459 
#> # A tibble: 5 × 12
#>   qname                             flag rname   pos  mapq cigar rnext pnext  tlen seq   qual  tags 
#>   <chr>                            <int> <chr> <int> <int> <chr> <chr> <int> <int> <chr> <chr> <chr>
#> 1 f0d7f308-026a-4e57-805f-73bfcc1…     4 *         0     0 *     *         0     0 TTAT… ")/&… ""   
#> 2 afbb49d7-926c-47b9-8e96-8a1c80b…     4 *         0     0 *     *         0     0 CTCT… "#$$… ""   
#> 3 82fbe615-2a10-4929-afd1-d09dc3c…     4 *         0     0 *     *         0     0 TTCA… "$%$… ""   
#> 4 314a7c4b-e91c-4704-95c8-b2ea9bb…     4 *         0     0 *     *         0     0 GGTA… "'''… ""   
#> 5 fab14603-7dda-4c01-a3d4-fd51718…     4 *         0     0 *     *         0     0 GGTG… "%%%… ""   
#> ... with 1715 more records

records(x)
#> # A tibble: 1,720 × 12
#>    qname                           flag rname    pos  mapq cigar rnext pnext  tlen seq   qual  tags 
#>    <chr>                          <int> <chr>  <int> <int> <chr> <chr> <int> <int> <chr> <chr> <chr>
#>  1 f0d7f308-026a-4e57-805f-73bfc…     4 *     0          0 *     *         0     0 TTAT… ")/&… ""   
#>  2 afbb49d7-926c-47b9-8e96-8a1c8…     4 *     0          0 *     *         0     0 CTCT… "#$$… ""   
#>  3 82fbe615-2a10-4929-afd1-d09dc…     4 *     0          0 *     *         0     0 TTCA… "$%$… ""   
#>  4 314a7c4b-e91c-4704-95c8-b2ea9…     4 *     0          0 *     *         0     0 GGTA… "'''… ""   
#>  5 fab14603-7dda-4c01-a3d4-fd517…     4 *     0          0 *     *         0     0 GGTG… "%%%… ""   
#>  6 60e737c0-8cb8-4ba6-8014-1f3a0…     0 chr12 1.23e8    60 213S… *         0     0 GGTG… "&%%… "NM:…
#>  7 1a953216-ddc0-4a24-9557-e3b7d…     0 chr2  3.29e7     1 343S… *         0     0 CAGT… "##\… "NM:…
#>  8 1a953216-ddc0-4a24-9557-e3b7d…   256 chr2  3.29e7     0 338H… *         0     0 GGGG… "%%$… "NM:…
#>  9 1a953216-ddc0-4a24-9557-e3b7d…   256 chr2  3.29e7     0 341H… *         0     0 GGGG… "%$%… "NM:…
#> 10 6998b387-7d93-4c5e-b323-28aed…    16 chr6  1.09e8    60 48S6… *         0     0 GGTT… "&(.… "NM:…
#> # ℹ 1,710 more rows
#> # ℹ Use `print(n = ...)` to see more rows
```

