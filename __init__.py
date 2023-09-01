#!/usr/bin/env python
# coding: utf-8

from ._finalreport import FinalReport


_FIELDS_ILLUMIN = ['SNP Name', 'Sample ID', 'Allele1 - AB', 'Allele2 - AB']
_RENAME_FIELDS = ['SNP_NAME', 'SAMPLE_ID', 'ALLELE1', 'ALLELE2']
_MAP_FIELDS = dict(zip(_FIELDS_ILLUMIN, _RENAME_FIELDS))

_ALLELE_CODE = {
	'AA': 0, 'AB': 1, 'BA': 1, 'BB': 2, '--': 5
}
