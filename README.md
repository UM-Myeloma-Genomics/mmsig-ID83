# mmsig Modification for ID83 Mutational Signature Analysis

## Overview

This repository provides a modification of the `mmsig` R package to perform mutational signature analysis using the **ID83** mutational signature. The ID83 mutational signatures include 83 types of insertion and deletion (indel) mutations, allowing for a more detailed analysis of mutational patterns.

## Installation

### 1. Install the `mmsig` Package

Ensure that you have the `mmsig` package installed in your R environment. You can install it from CRAN:

```R
# install.packages("devtools")
devtools::install_github("evenrus/mmsig")
```

### 2. Import the Modified `mmsig_ID83.R` Script

Clone or download this repository, and source the `mmsig_ID83.R` script into your R environment:

```R
source("path/to/mmsig_ID83.R")
```

## Usage

### Input Data Format

- **ID83 Matrix Format**: Your input data must be in the ID83 matrix format, which includes counts of 83 types of indel mutations.
- **Generating the Matrix**: For easy matrix creation, you can use the [SigProfiler Matrix Generator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator) tool. Use the `.all` output file, which contains the necessary indel counts.

### Example

```R
# Load required library
library(mmsig)

# Source the modified mmsig script
source("path/to/mmsig_ID83.R")

# Load your ID83 formatted data
id83_data <- read.table("path/to/your_data.all", header = TRUE, sep = "\t", row.names = 1)

# Perform mutational signature analysis
results <- mmsig_ID83(id83_data)

# View the results
print(results)
```

## Notes

- **Data Preparation**: Ensure your input data is correctly formatted according to the ID83 standard.
- **Script Modification**: The `mmsig_ID83.R` script contains modifications to handle ID83 data specifically.

