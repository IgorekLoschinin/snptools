# snptools

**Snptools** is a tool for SNP (Single Nucleotide Polymorphism) data processing, 
parentage calculation and call rate estimation.

## Introduction

SNP (Single Nucleotide Polymorphism) represent genetic variations, that can 
be used to analyze genetic data. SNPTools provides a set of tools for working 
with SNP data, including the following capabilities:

- SNP data processing - FinalReport.
- Parentage Verification and Parentage Discovery Based on SNP Genotypes (ICAR). 
- Call rate estimation (percentage of missing data).
- Processing and preparation of data in plink formats.

## Installation

To install SNPTools, follow the steps below:

1. Clone the repository into your project directory:
   ```
   git clone https://github.com/yourusername/snpTools.git
   ```
2. Set dependencies:
   ```
   pip install -r requirements.txt
   ```
3. Use SNPTools:
   ```
   import snptools
   ```
   
## Usage
Snptools provides commands for a variety of operations. Here are examples of 
usage:

#### SNP data processing:
```
from snptools.finalreport import FinalReport
```

#### Computation of parentage:
```
from snptools.parentage import Discovery, Verification
```

#### Preparation format files:
```
from snptools.format import (
   Snp, make_fam, make_ped, make_lgen, make_map
)
```

#### Stat:
```
from snptools.statistics import (
   hwe, hwe_test, call_rate, allele_freq, minor_allele_freq
)
```

## Documentation
Detailed documentation on how to use SNPTools is available see the docs.

## License
This project is licensed under the GNU General Public License - see the 
LICENSE file for details.
