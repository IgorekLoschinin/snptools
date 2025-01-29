#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

__all__ = [
    "make_map", "make_ped", "make_fam", "make_lgen"
]

import re
import pandas as pd


def make_map(manifest: pd.DataFrame) -> pd.DataFrame:
    """ PLINK text fileset variant information file
    https://www.cog-genomics.org/plink/1.9/formats#map

    A text file with no header line, and one line per variant with the
    following 3-4 fields:

        1. Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
        2. Variant identifier.
        3. Position in morgans or centimorgans (optional; also safe to use dummy value of '0').
        4. Base-pair coordinate.

    All lines must have the same number of columns (so either no lines contain
    the morgans/centimorgans column, or all of them do).

    :param manifest: The file that is taken on the Illumina website with full information about the chip https://support.illumina.com/downloads/bovinesnp50-v3-0-product-files.html.

    :return: Return data in formate .map.
    """

    fields = ['Chr', 'Name', 'MapInfo']

    if all([
        True
        if item not in manifest.columns
        else False
        for item in fields
    ]):
        raise KeyError("Manifest has no data to build map format!")

    # Rearrange the columns and replace the names of the sex and mitochondrial
    # chromosomes
    permute_cols = manifest[fields].\
        sort_values(by='Name').\
        replace({'X': 30, 'Y': 31, 'MT': 33}).\
        dropna(axis=0)

    # Insert distances in centimorganides
    permute_cols.insert(2, 'morgans', [0] * len(manifest))

    return permute_cols


def make_ped(
        data: pd.DataFrame,
        sid_col: str,
        snp_col: str,
        fid_col: str = None,
        father_col: str = None,
        mother_col: str = None,
        sex_col: str = None,
) -> pd.DataFrame | None:
    """ Original standard text format for sample pedigree information and
    genotype calls. Normally must be accompanied by a .map file.
    https://www.cog-genomics.org/plink/1.9/formats#ped

    The PED file has 6 fixed columns at the beginning followed by the SNP
    information. The columns should be separated by a whitespace or a tab. The
    first six columns hold the following information:

        1. Family ID (if unknown use the same id as for the sample id in column two)
        2. Sample ID
        3. Paternal ID (if unknown use 0)
        4. Maternal ID (if unknown use 0)
        5. Sex (1=male; 2=female; 0=unknown)
        6. Affection (0=unknown; 1=unaffected; 2=affected)
        7. Genotypes (space or tab separated, 2 for each marker. 0/-9=missing)

    Here is a brief example of a genotype PED file containing 5 samples
    with 10 homozygous SNPs:

    4304 4304 0 0 0 0 C C C C G G G G G G C C G G C C T T T T
    6925 6925 0 0 0 0 C C C C T T G G A A C C G G C C T T T T
    7319 7319 0 0 0 0 C C C C G G G G G G C C G G C C T T T T
    6963 6963 0 0 0 0 A A C C T T G G A A C C G G C C T T T T
    6968 6968 0 0 0 0 C C C C G G G G G G G G G G C C T T T T

    :param data: Snp data that contain full or partial information on the animal.
    :param sid_col: Sample ID. Column name in data.
    :param snp_col: Snp column name in data.
    :param fid_col: Family ID column name in data (if unknown use the same id as for the sample id in column two).
    :param father_col: Paternal ID column name in data (if unknown use 0).
    :param mother_col: Maternal ID column name in data (if unknown use 0).
    :param sex_col: Sex column name in data (if unknown use 0).
    :return: Returns an array of data in ped format to work with the plink program.
    """

    _fields = ["fid", "sid", "father", "mother", "sex", "not_used", "snp"]
    _f_dtype = dict(zip(_fields, (str for _ in range(len(_fields)))))

    _ped = pd.DataFrame(columns=_fields)

    if sid_col not in data.columns or snp_col not in data.columns:
        raise KeyError(f"Data has not in name columns!")

    # Checked Sample ID on underscope - '_'
    _ped["sid"] = data[sid_col].astype(str)
    if _ped["sid"].apply(_check_underscore).any():
        raise Exception(
            "Replace in 'Sample ID' columns '_' on another a simbols"
        )

    # Checked Family ID on underscope - '_'
    if fid_col is not None:
        if fid_col not in data.columns:
            raise KeyError(f"Data has not in name columns {fid_col}!")

        if (data[fid_col].dtype.hasobject and
                data[fid_col].apply(_check_underscore).any()):
            raise Exception(
                "Replace in 'Family ID' columns '_' on another a simbols"
            )

        _ped["fid"] = data[fid_col]

    else:
        _ped["fid"] = data[sid_col].astype(str)

    _ped["father"] = data[father_col] if father_col is not None else 0
    _ped["mother"] = data[mother_col] if mother_col is not None else 0
    _ped["sex"] = data[sex_col] if sex_col is not None else 0
    _ped["not_used"] = 0
    _ped["snp"] = data[snp_col]

    return _ped[_fields].astype(_f_dtype)


def make_fam(
        data: pd.DataFrame,
        sid_col: str,
        fid_col: str = None,
        father_col: str = None,
        mother_col: str = None,
        sex_col: str = None,
        sex_val: int = 0,
        pheno_col: str = None,
        pheno_val: int = -9
) -> pd.DataFrame | None:
    """ PLINK sample information file
    https://www.cog-genomics.org/plink/1.9/formats#fam

    A text file with no header line, and one line per sample with the following
    six fields:

        1. Family ID ('FID')
        2. Within-family ID ('IID'; cannot be '0')
        3. Within-family ID of father ('0' if father isn't in dataset)
        4. Within-family ID of mother ('0' if mother isn't in dataset)
        5. Sex code ('1' = male, '2' = female, '0' = unknown)
        6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

    :param data: Snp data that contain full or partial information on the animal
    :param fid_col: Family ID, default value "1". Must not contain underline - "_"
    :param sid_col: Within-family ID ('IID'; cannot be '0'). Must not contain underline - "_"
    :param father_col: Within-family ID of father ('0' if father isn't in dataset)
    :param mother_col: Within-family ID of mother ('0' if mother isn't in dataset)
    :param sex_col: Sex column name in data
    :param sex_val: Sex code ('1' = male, '2' = female, '0' = unknown)
    :param pheno_col: Pheno column name in data
    :param pheno_val: Phenotype value ('1' = control, '2' = case,'-9'/'0'/non-numeric = missing data if case/control)
    :return: Return data in formate .fam
    """

    _fields = ['fid', 'sid', 'father', 'mother', 'sex', 'pheno']
    _f_dtype = dict(zip(_fields, (str for _ in range(len(_fields)))))

    _fam = pd.DataFrame(columns=_fields)

    if sid_col not in data.columns:
        raise KeyError(f"Data has not in name columns {sid_col}!")

    # Checked Sample ID on underscope - '_'
    _fam["sid"] = data[sid_col].astype(str)
    if _fam["sid"].apply(_check_underscore).any():
        raise Exception(
            "Replace in 'Sample ID' columns '_' on another a simbols"
        )

    # Checked Family ID on underscope - '_'
    if fid_col is not None:
        if fid_col not in data.columns:
            raise KeyError(f"Data has not in name columns {fid_col}!")

        if (data[fid_col].dtype.hasobject and
                data[fid_col].apply(_check_underscore).any()):
            raise Exception(
                "Replace in 'Family ID' columns '_' on another a simbols"
            )

        _fam["fid"] = data[fid_col]

    else:
        _fam["fid"] = 1

    _fam["father"] = data[father_col] if father_col is not None else 0
    _fam["mother"] = data[mother_col] if mother_col is not None else 0
    _fam["sex"] = data[sex_col] if sex_col is not None else sex_val
    _fam['pheno'] = data[pheno_col] if pheno_col is not None else pheno_val

    return _fam[_fields].astype(_f_dtype)


def make_lgen(
        data: pd.DataFrame,
        sid_col: str,
        snp_name: str,
        alleles: list[str],
        fid_col: str = None
) -> pd.DataFrame | None:
    """ PLINK long-format genotype file
    https://www.cog-genomics.org/plink/1.9/formats#lgen

    A text file with no header line, and one line per genotype call (or
    just not-homozygous-major calls if 'lgen-ref' was invoked) usually with
    the following five fields:

        1. Family ID
        2. Within-family ID
        3. Variant identifier
        4. Allele call 1 ('0' for missing)
        5. Allele call 2

    There are several variations which are also handled by PLINK; see the
    original discussion for details.

    :param data: Data the after parsing FinalReport.txt
    :param sid_col:
    :param snp_name:
    :param fid_col: Family ID, default value "1"
    :param alleles:
    :return: Return data in formate .lgen
    """
    _fields = ['fid', 'sid', 'snp_name', 'allele1', 'allele2']
    _f_dtype = dict(zip(_fields, (str for _ in range(len(_fields)))))

    _lgen = pd.DataFrame(columns=_fields)

    try:
        # Checked Sample ID on underscope - '_'
        _lgen["sid"] = data[sid_col].astype(str)
        if _lgen["sid"].apply(_check_underscore).any():
            raise Exception(
                "Replace in 'Sample ID' columns '_' on another a simbols"
            )

        # Checked Family ID on underscope - '_'
        if fid_col is not None:
            if (data[fid_col].dtype.hasobject and
                    data[fid_col].apply(_check_underscore).any()):
                raise Exception(
                    "Replace in 'Family ID' columns '_' on another a simbols"
                )

            _lgen["fid"] = data[fid_col]

        else:
            _lgen["fid"] = 1

        _lgen["snp_name"] = data[snp_name]
        _lgen[["allele1", "allele2"]] = data[alleles].replace({'-': 0})

    except Exception as e:
        raise e

    return _lgen[_fields].astype(_f_dtype)


def _check_underscore(value: str) -> bool:
    """ Checking for underscore in a string

    :param value: String for checked
    :return: Return True if there is an underline in the string, False if not.
    """
    _under_l = re.compile(r"_")

    if _under_l.findall(value):
        return True

    return False
