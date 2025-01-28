Usage
=====

Snptools provides commands for a variety of operations. Here are examples of
usage.

**SNP data processing**::

    from snplib.finalreport import FinalReport


**Computation of parentage**::

    from snplib.parentage import Discovery, Verification

**Preparation format files**::

    from snplib.format import (
        Snp, make_fam, make_ped, make_lgen, make_map
    )

**Stat**::

    from snplib.statistics import (
       hwe, hwe_test, call_rate, allele_freq, minor_allele_freq
    )
