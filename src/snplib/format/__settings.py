#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

FIELDS_ILLUMIN = ['SNP Name', 'Sample ID', 'Allele1 - AB', 'Allele2 - AB']
RENAME_FIELDS = ['SNP_NAME', 'SAMPLE_ID', 'ALLELE1', 'ALLELE2']
MAP_FIELDS = dict(zip(FIELDS_ILLUMIN, RENAME_FIELDS))
