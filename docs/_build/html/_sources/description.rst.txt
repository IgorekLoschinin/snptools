Introduction
============
**Snptools** is a tool for SNP (Single Nucleotide Polymorphism) data processing,
parentage calculation and call rate estimation.

SNP (Single Nucleotide Polymorphism) represent genetic variations, that can 
be used to analyze genetic data. SNPTools provides a set of tools for working 
with SNP data, including the following capabilities:

- Processing the Finalreport.txt file in Dataframe.
- Preparation of data in plink and snp(blupf90) formats.
- Parentage Verification and Parentage Discovery Based on SNP Genotypes (ICAR). 
- Call rate estimation (percentage of missing data).
- Allele frequency, minor allele frequency (MAF).
- Hardy-Weinberg equilibrium.

Installation
============

To install SNPTools, follow the steps below:

Clone the repository into your project directory::

   git clone https://github.com/yourusername/snpTools.git

Set dependencies::

   python -m pip install -r requirements.txt

Use SNPTools::

   import snptools

Usage
=====

Snptools provides commands for a variety of operations. Here are examples of
usage:

SNP data processing::

    from snptools.finalreport import FinalReport


Computation of parentage::

    from snptools.parentage import Discovery, Verification

Preparation format files::

    from snptools.format import (
        Snp, make_fam, make_ped, make_lgen, make_map
    )

Stat::

    from snptools.statistics import (
       hwe, hwe_test, call_rate, allele_freq, minor_allele_freq
    )
