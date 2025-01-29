# snptools
<p align="center">
  <img width="150" height="150" src="./iconlib.png">
</p>

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
You can install snptools via pip from [PyPI](https://pypi.org/project/snplib/):
```
pip install snplib
```

## Usage
Snptools provides commands for a variety of operations. Here are examples of 
usage:

#### SNP data processing:
```
from snplib import FinalReport
```

#### Computation of parentage:
```
from snplib import Discovery, Verification
```

#### Preparation format files:
```
from snplib import (
   Snp, make_fam, make_ped, make_lgen, make_map
)
```

#### Stat:
```
from snplib import (
   hwe, hwe_test, call_rate, allele_freq, minor_allele_freq
)
```

## Documentation
Detailed documentation on how to use SNPTools is available see the [docs](https://snptools.readthedocs.io/en/latest/).

## License
This project is licensed under the GNU General Public License - see the 
LICENSE file for details.
