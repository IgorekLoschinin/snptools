Code Examples
=============

Head file FinalReport
---------------------
::

    [Header]
    GSGT Version	2.0.4
    Processing Date	1/16/2023 9:19 AM
    Content		GGP_HDv3_C.bpm
    Num SNPs	139376
    Total SNPs	140668
    Num Samples	25
    Total Samples	26
    File 	1 of 25
    [Data]
    SNP Name	Sample ID	Allele1 - Forward	Allele2 - Forward	Allele1 - Top	Allele2 - Top	Allele1 - AB	Allele2 - AB	GC Score	X	Y
    ARS-BFGL-BAC-10172	HO840M003135245650	G	G	G	G	B	B	0.9420	0.069	0.801
    ARS-BFGL-BAC-1020	HO840M003135245650	G	G	G	G	B	B	0.9489	0.033	0.700
    ARS-BFGL-BAC-10245	HO840M003135245650	C	C	G	G	B	B	0.7277	0.152	1.504
    ARS-BFGL-BAC-10345	HO840M003135245650	A	C	A	C	A	B	0.9411	0.598	0.572
    ARS-BFGL-BAC-10375	HO840M003135245650	A	G	A	G	A	B	0.9348	0.430	0.494
    ...

Processing the Finalreport.txt
------------------------------

.. code-block:: python

        from snplib.finalreport import FinalReport

        obj_report = FinalReport()
        obj_report.handle("path/to/finalreport.txt")

        header = obj_report.header
        data = obj_report.snp_data

Output Dataframe::

    SNP Name	Sample ID	Allele1 - Forward	Allele2 - Forward	Allele1 - Top	Allele2 - Top	Allele1 - AB	Allele2 - AB	GC Score	X	Y
    ARS-BFGL-BAC-10172	HO840M003135245650	G	G	G	G	B	B	0.9420	0.069	0.801
    ARS-BFGL-BAC-1020	HO840M003135245650	G	G	G	G	B	B	0.9489	0.033	0.700
    ARS-BFGL-BAC-10245	HO840M003135245650	C	C	G	G	B	B	0.7277	0.152	1.504
    ARS-BFGL-BAC-10345	HO840M003135245650	A	C	A	C	A	B	0.9411	0.598	0.572
    ARS-BFGL-BAC-10375	HO840M003135245650	A	G	A	G	A	B	0.9348	0.430	0.494


Select alleles - AB or Top or Forward

.. code-block:: python

        alleles_ab = FinalReport(allele="AB")
        alleles_ab.handle("path/to/finalreport.txt")
        data_ab = alleles_ab.snp_data

        alleles_top = FinalReport(allele="Top")
        alleles_top.handle("path/to/finalreport.txt")
        data_top = alleles_top.snp_data

        alleles_forward = FinalReport(allele="Forward")
        alleles_forward.handle("path/to/finalreport.txt")
        data_forward = alleles_forward.snp_data

Output::

        SNP Name    Sample ID    Allele1 - AB  Allele2 - AB     GC Score    X Y
        ARS-BFGL-BAC-10172 HO840M003135245650 B B 0.9420 0.069 0.801
        ARS-BFGL-BAC-1020 HO840M003135245650 B B 0.9489 0.033 0.700
        ARS-BFGL-BAC-10245 HO840M003135245650 B B 0.7277 0.152 1.504
        ARS-BFGL-BAC-10345 HO840M003135245650 A B 0.9411 0.598 0.572
        ARS-BFGL-BAC-10375 HO840M003135245650 A B 0.9348 0.430 0.494f

        ...

To handle large files, use `usecols` and `dtype`. This reduces memory
consumption and speeds up processing.

.. note::
    `usecols` is used when `allele` is **None**.

.. code-block:: python

        alleles_ab = FinalReport(
            usecols=['SNP Name', 'Sample ID', 'Allele1 - AB', 'Allele2 - AB'],
            dtype={'SNP Name': 'category'}
        )
        alleles_ab.handle("path/to/finalreport.txt")
        data_ab = alleles_ab.snp_data

Output::

        SNP Name    Sample ID    Allele1 - AB  Allele2 - AB
        ARS-BFGL-BAC-10172 HO840M003135245650 B B
        ARS-BFGL-BAC-1020 HO840M003135245650 B B
        ARS-BFGL-BAC-10245 HO840M003135245650 B B
        ARS-BFGL-BAC-10345 HO840M003135245650 A B
        ARS-BFGL-BAC-10375 HO840M003135245650 A B
        ...

Preparation SNP files
---------------------

After processing the raw data, FinalReports.txt, for further analysis
several steps of SNP (Single Nucleotide Polymorphism) file preparation are
necessary.

Data formatting
---------------

The received data often requires formatting to bring it to a standardized form.
The proposed module includes data formatting for the programs blupf90 and
plink - GBLUP, ssGBLUP, GWAS.

blupf90 format
______________
The input data for obtaining the ``snp.txt`` file used for the genomic
blupf90 evaluation is the data file - processed file ``finalreport.txt``.
The processed file can be seen in the item above - Finalreport.txt processing:

Content input *file.txt*::

        SNP Name    Sample ID    Allele1 - AB  Allele2 - AB     GC Score    X Y
        ARS-BFGL-BAC-10172 HO840M003135245650 B B 0.9420 0.069 0.801
        ARS-BFGL-BAC-1020 HO840M003135245650 B B 0.9489 0.033 0.700
        ARS-BFGL-BAC-10245 HO840M003135245650 B B 0.7277 0.152 1.504
        ARS-BFGL-BAC-10345 HO840M003135245650 A B 0.9411 0.598 0.572
        ARS-BFGL-BAC-10375 HO840M003135245650 A B 0.9348 0.430 0.494f

        ...

**uga**

.. code-block:: python

    import pandas as pd
    from snplib.format import Snp

    data_finalreport = pd.read_csv("path_to_file/file.txt", sep="\t")

    obj = Snp(fmt="uga")
    obj_snp.process(data_finalreport)
    obj_snp.to_file("./snp.txt")

Data after snp processing in ``uga`` (blupf90) format - obj_snp.data::

      SAMPLE_ID                SNP
    0     14814  02011015010000500
    1     14815  01110152120222512

Default result - this is what the data looks like if ``fmt=None``::

                    SNP_NAME SAMPLE_ID SNP
    0               ABCA12     14814   0
    1   ARS-BFGL-BAC-13031     14814   2
    2   ARS-BFGL-BAC-13039     14814   0
    3   ARS-BFGL-BAC-13049     14814   1
                    ...
    17              ABCA12     14815   0
    18  ARS-BFGL-BAC-13031     14815   1
    19  ARS-BFGL-BAC-13039     14815   1
    20  ARS-BFGL-BAC-13049     14815   1
                    ...

plink format
____________

This page describes specialized PLINK input and output file formats which are
identifiable by file extension. https://www.cog-genomics.org/plink/1.9/formats
Common fomrats for performing GWAS analysis - ``ped``, ``map``, ``fam``, ``lgen``....

**map** - https://www.cog-genomics.org/plink/1.9/formats#map

To get the ``.map`` file, first you need to download the *manifest file* for the chip
you are using chip.

.. note::
    *file_bovinesnp50.csv* - The file that is taken on the Illumina website with full
    information about the chip https://support.illumina.com/downloads/bovinesnp50-v3-0-product-files.html

Since the make_map function accepts **pd.DataFrame**, the *manifest file* processing is performed
independently.

Input data for make_map::

                                           IlmnID  ... BeadSetID
    0       BovineHD0100037694-128_T_F_2278925834  ...      1241
    1   BovineHD0100037699_dup-128_T_F_2327674593  ...      1241
    2   BovineHD0100037703_dup-128_B_R_2327674602  ...      1241
    3   BovineHD0100037704_dup-128_T_F_2327674603  ...      1241
    4   BovineHD0100037710_dup-128_T_F_2327674613  ...      1241
    5   BovineHD0100037712_dup-128_B_R_2327674618  ...      1241
    6       BovineHD0100037716-128_T_F_2255347065  ...      1241
    7       BovineHD0100037719-128_T_F_2278926219  ...      1241
    8       BovineHD0100037720-128_B_R_2255342455  ...      1241
    9   BovineHD0100037722_dup-128_B_R_2327674634  ...      1241


.. note::
    The original file, for example, **BovineSNP50_v3_A1.csv** looks like this::

        Illumina, Inc.,,,,,,,,,,,,,,,,,
        [Heading],,,,,,,,,,,,,,,,,,
        Descriptor File Name,BovineSNP50_v3_A1.bpm,,,,,,,,,,,,,,,,,
        Assay Format,Infinium HTS,,,,,,,,,,,,,,,,,
        Date Manufactured,1/14/2016,,,,,,,,,,,,,,,,,
        Loci Count ,53218,,,,,,,,,,,,,,,,,
        [Assay],,,,,,,,,,,,,,,,,,
        IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID
        ABCA12_r2-1_T_F_2277749139,ABCA12,TOP,[A/G],0059616496,CTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC,,,0,2,103548215,diploid,Bos taurus,UMD3.1,1,TOP,ACTCTGGTGGATGGTTCATAATCTGCTAAGATGAATAAGTTACTGGGGAAACTGGTGCATTTATTTTAAATATAAATTATATAGTCTGTAAGATATAAAGACTGCCTAATTTATTTGAACACCATACTGATCTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC[A/G]CTCCCTTCCAGCTTACCTCAACAGCCTGAATAATTTCCTCCTGCGAGTTAACATGTCAAAATATGATGCTGCCCGACATGGTAAAGTTATTTACATAGGAGCTCCTTGTATTGAAACTCTTGCTACTCTCCATGTGAAAATATACATTAGACCCCATTTTCCTCCCTGTGGCAGCTAT,ACTCTGGTGGATGGTTCATAATCTGCTAAGATGAATAAGTTACTGGGGAAACTGGTGCATTTATTTTAAATATAAATTATATAGTCTGTAAGATATAAAGACTGCCTAATTTATTTGAACACCATACTGATCTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC[A/G]CTCCCTTCCAGCTTACCTCAACAGCCTGAATAATTTCCTCCTGCGAGTTAACATGTCAAAATATGATGCTGCCCGACATGGTAAAGTTATTTACATAGGAGCTCCTTGTATTGAAACTCTTGCTACTCTCCATGTGAAAATATACATTAGACCCCATTTTCCTCCCTGTGGCAGCTAT,1241
        APAF1_dup-1_B_F_2327661418,APAF1,BOT,[T/C],0041654401,ATATTGTGCAACTGGGCCTCTGTGAACTGGAAACTTCAGAGGTTTATCGG,,,0,5,63150400,diploid,Bos taurus,UMD3.1,1,BOT,CCATTTCCTAATATTGTGCAACTGGGCCTCTGTGAACTGGAAACTTCAGAGGTTTATCGG[T/C]AAGCTAAGCTGCAGGCCAAGCAGGAGGTCGATAACGGAATGCTTTACCTGGAGTGGGTGT,ACACCCACTCCAGGTAAAGCATTCCGTTATCGACCTCCTGCTTGGCCTGCAGCTTAGCTT[A/G]CCGATAAACCTCTGAAGTTTCCAGTTCACAGAGGCCCAGTTGCACAATATTAGGAAATGG,1241
        ARS-BFGL-BAC-10172_dup-0_T_F_2328966397,ARS-BFGL-BAC-10172,TOP,[A/G],0072620471,GGTCCCCAAAGTATGTGGTAGCACTTACTTATGTAAGTCATCACTCAAGT,,,3,14,6371334,diploid,Bos taurus,UM3,0,TOP,CTCAGAAGTTGGTCCCCAAAGTATGTGGTAGCACTTACTTATGTAAGTCATCACTCAAGT[A/G]ATCCAGAATATTCTTTTAGTAATATTTTTGTTAATATTGAAATTTTTAAAACAATTGAAA,CTCAGAAGTTGGTCCCCAAAGTATGTGGTAGCACTTACTTATGTAAGTCATCACTCAAGT[A/G]ATCCAGAATATTCTTTTAGTAATATTTTTGTTAATATTGAAATTTTTAAAACAATTGAAA,1241
        .
        .
        .
        UA-IFASA-9812_dup-0_B_F_2329051536,UA-IFASA-9812,BOT,[T/C],0031677304,ACCTCCATAGCTGATAGGAATGGTCTCAACTTGCAGCCCCATTATACTAA,,,3,29,48012818,diploid,Bos taurus,UM3,0,BOT,GTAAAAACAAACCTCCATAGCTGATAGGAATGGTCTCAACTTGCAGCCCCATTATACTAA[T/C]GATGATCTGAAGTTTCTCAAGCACGCAGAGAAACGTAAGAGAAACGTTCCAGCAAAGGGA,TCCCTTTGCTGGAACGTTTCTCTTACGTTTCTCTGCGTGCTTGAGAAACTTCAGATCATC[A/G]TTAGTATAATGGGGCTGCAAGTTGAGACCATTCCTATCAGCTATGGAGGTTTGTTTTTAC,1241
        UA-IFASA-9813_dup-0_B_F_2329051538,UA-IFASA-9813,BOT,[T/C],0011661313,ACCTTTGCACTCGCTAACGGTTCAGCATTAATCAGACTTCCTCAGGAATT,,,3,19,32508700,diploid,Bos taurus,UM3,0,BOT,AATAAAACCAACCTTTGCACTCGCTAACGGTTCAGCATTAATCAGACTTCCTCAGGAATT[T/C]AGGGGTCAATTCCCCCATGTCTAAAATTGAACCTCAACGTCCTTTCTGTTTTCAAAACTC,GAGTTTTGAAAACAGAAAGGACGTTGAGGTTCAATTTTAGACATGGGGGAATTGACCCCT[A/G]AATTCCTGAGGAAGTCTGATTAATGCTGAACCGTTAGCGAGTGCAAAGGTTGGTTTTATT,1241
        UMPS_dup-1_T_R_2327737250,UMPS,TOP,[A/G],0073777348,TAACTGAACTCCTGGAGTCAAGTGAAGAAATTCTGGTTTCATGCTTACTC,,,0,1,69756880,diploid,Bos taurus,UMD3.1,1,BOT,TCATCTGTTGATTACATTCCATTCAGGTGCAAATGGCTGAAGAACATTCTGAATTTGTGATTGGTTTTATTTCTGGCTCC[T/C]GAGTAAGCATGAAACCAGAATTTCTTCACTTGACTCCAGGAGTTCAGTTAGAAGCAGGAGGTAAGCCTATTGATTGGTAA,TTACCAATCAATAGGCTTACCTCCTGCTTCTAACTGAACTCCTGGAGTCAAGTGAAGAAATTCTGGTTTCATGCTTACTC[A/G]GGAGCCAGAAATAAAACCAATCACAAATTCAGAATGTTCTTCAGCCATTTGCACCTGAATGGAATGTAATCAACAGATGA,1241
        [Controls],,,,,,,,,,,,,,,,,,
        0027630314,Staining,Red,DNP (High),,,,,,,,,,,,,,,
        0029619375,Staining,Purple,DNP (Bgnd),,,,,,,,,,,,,,,
        0041666334,Staining,Green,Biotin (High),,,,,,,,,,,,,,,
        0034648333,Staining,Blue,Biotin (Bgnd),,,,,,,,,,,,,,,
        0017616306,Extension,Red,Extension (A),,,,,,,,,,,,,,,
        0014607337,Extension,Purple,Extension (T),,,,,,,,,,,,,,,

    Therefore, for direct reading via **pd.read_csv()** it is necessary to
    preprocess the file - delete extra lines::

        Illumina, Inc.,,,,,,,,,,,,,,,,,
        [Heading],,,,,,,,,,,,,,,,,,
        Descriptor File Name,BovineSNP50_v3_A1.bpm,,,,,,,,,,,,,,,,,
        Assay Format,Infinium HTS,,,,,,,,,,,,,,,,,
        Date Manufactured,1/14/2016,,,,,,,,,,,,,,,,,
        Loci Count ,53218,,,,,,,,,,,,,,,,,
        [Assay],,,,,,,,,,,,,,,,,,

        and

        [Controls],,,,,,,,,,,,,,,,,,
        0027630314,Staining,Red,DNP (High),,,,,,,,,,,,,,,
        0029619375,Staining,Purple,DNP (Bgnd),,,,,,,,,,,,,,,
        0041666334,Staining,Green,Biotin (High),,,,,,,,,,,,,,,
        0034648333,Staining,Blue,Biotin (Bgnd),,,,,,,,,,,,,,,
        0017616306,Extension,Red,Extension (A),,,,,,,,,,,,,,,
        0014607337,Extension,Purple,Extension (T),,,,,,,,,,,,,,,

    The file should end up looking like this::

        IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID
        ABCA12_r2-1_T_F_2277749139,ABCA12,TOP,[A/G],0059616496,CTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC,,,0,2,103548215,diploid,Bos taurus,UMD3.1,1,TOP,ACTCTGGTGGATGGTTCATAATCTGCTAAGATGAATAAGTTACTGGGGAAACTGGTGCATTTATTTTAAATATAAATTATATAGTCTGTAAGATATAAAGACTGCCTAATTTATTTGAACACCATACTGATCTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC[A/G]CTCCCTTCCAGCTTACCTCAACAGCCTGAATAATTTCCTCCTGCGAGTTAACATGTCAAAATATGATGCTGCCCGACATGGTAAAGTTATTTACATAGGAGCTCCTTGTATTGAAACTCTTGCTACTCTCCATGTGAAAATATACATTAGACCCCATTTTCCTCCCTGTGGCAGCTAT,ACTCTGGTGGATGGTTCATAATCTGCTAAGATGAATAAGTTACTGGGGAAACTGGTGCATTTATTTTAAATATAAATTATATAGTCTGTAAGATATAAAGACTGCCTAATTTATTTGAACACCATACTGATCTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC[A/G]CTCCCTTCCAGCTTACCTCAACAGCCTGAATAATTTCCTCCTGCGAGTTAACATGTCAAAATATGATGCTGCCCGACATGGTAAAGTTATTTACATAGGAGCTCCTTGTATTGAAACTCTTGCTACTCTCCATGTGAAAATATACATTAGACCCCATTTTCCTCCCTGTGGCAGCTAT,1241
        APAF1_dup-1_B_F_2327661418,APAF1,BOT,[T/C],0041654401,ATATTGTGCAACTGGGCCTCTGTGAACTGGAAACTTCAGAGGTTTATCGG,,,0,5,63150400,diploid,Bos taurus,UMD3.1,1,BOT,CCATTTCCTAATATTGTGCAACTGGGCCTCTGTGAACTGGAAACTTCAGAGGTTTATCGG[T/C]AAGCTAAGCTGCAGGCCAAGCAGGAGGTCGATAACGGAATGCTTTACCTGGAGTGGGTGT,ACACCCACTCCAGGTAAAGCATTCCGTTATCGACCTCCTGCTTGGCCTGCAGCTTAGCTT[A/G]CCGATAAACCTCTGAAGTTTCCAGTTCACAGAGGCCCAGTTGCACAATATTAGGAAATGG,1241
        ARS-BFGL-BAC-10172_dup-0_T_F_2328966397,ARS-BFGL-BAC-10172,TOP,[A/G],0072620471,GGTCCCCAAAGTATGTGGTAGCACTTACTTATGTAAGTCATCACTCAAGT,,,3,14,6371334,diploid,Bos taurus,UM3,0,TOP,CTCAGAAGTTGGTCCCCAAAGTATGTGGTAGCACTTACTTATGTAAGTCATCACTCAAGT[A/G]ATCCAGAATATTCTTTTAGTAATATTTTTGTTAATATTGAAATTTTTAAAACAATTGAAA,CTCAGAAGTTGGTCCCCAAAGTATGTGGTAGCACTTACTTATGTAAGTCATCACTCAAGT[A/G]ATCCAGAATATTCTTTTAGTAATATTTTTGTTAATATTGAAATTTTTAAAACAATTGAAA,1241
        .
        .
        .
        UA-IFASA-9812_dup-0_B_F_2329051536,UA-IFASA-9812,BOT,[T/C],0031677304,ACCTCCATAGCTGATAGGAATGGTCTCAACTTGCAGCCCCATTATACTAA,,,3,29,48012818,diploid,Bos taurus,UM3,0,BOT,GTAAAAACAAACCTCCATAGCTGATAGGAATGGTCTCAACTTGCAGCCCCATTATACTAA[T/C]GATGATCTGAAGTTTCTCAAGCACGCAGAGAAACGTAAGAGAAACGTTCCAGCAAAGGGA,TCCCTTTGCTGGAACGTTTCTCTTACGTTTCTCTGCGTGCTTGAGAAACTTCAGATCATC[A/G]TTAGTATAATGGGGCTGCAAGTTGAGACCATTCCTATCAGCTATGGAGGTTTGTTTTTAC,1241
        UA-IFASA-9813_dup-0_B_F_2329051538,UA-IFASA-9813,BOT,[T/C],0011661313,ACCTTTGCACTCGCTAACGGTTCAGCATTAATCAGACTTCCTCAGGAATT,,,3,19,32508700,diploid,Bos taurus,UM3,0,BOT,AATAAAACCAACCTTTGCACTCGCTAACGGTTCAGCATTAATCAGACTTCCTCAGGAATT[T/C]AGGGGTCAATTCCCCCATGTCTAAAATTGAACCTCAACGTCCTTTCTGTTTTCAAAACTC,GAGTTTTGAAAACAGAAAGGACGTTGAGGTTCAATTTTAGACATGGGGGAATTGACCCCT[A/G]AATTCCTGAGGAAGTCTGATTAATGCTGAACCGTTAGCGAGTGCAAAGGTTGGTTTTATT,1241
        UMPS_dup-1_T_R_2327737250,UMPS,TOP,[A/G],0073777348,TAACTGAACTCCTGGAGTCAAGTGAAGAAATTCTGGTTTCATGCTTACTC,,,0,1,69756880,diploid,Bos taurus,UMD3.1,1,BOT,TCATCTGTTGATTACATTCCATTCAGGTGCAAATGGCTGAAGAACATTCTGAATTTGTGATTGGTTTTATTTCTGGCTCC[T/C]GAGTAAGCATGAAACCAGAATTTCTTCACTTGACTCCAGGAGTTCAGTTAGAAGCAGGAGGTAAGCCTATTGATTGGTAA,TTACCAATCAATAGGCTTACCTCCTGCTTCTAACTGAACTCCTGGAGTCAAGTGAAGAAATTCTGGTTTCATGCTTACTC[A/G]GGAGCCAGAAATAAAACCAATCACAAATTCAGAATGTTCTTCAGCCATTTGCACCTGAATGGAATGTAATCAACAGATGA,1241


.. code-block:: python

    import pandas as pd
    from snplib.format import make_map

    input_data = pd.read_csv("./file_bovinesnp50.csv")
    data_map = make_map(input_data)

Output data view::

        Chr                Name  morgans  MapInfo
         0  BovineHD0100037694        0        0
         0  BovineHD0100037699        0        0
         0  BovineHD0100037703        0        0
         0  BovineHD0100037704        0        0

**ped** - https://www.cog-genomics.org/plink/1.9/formats#ped

.. code-block:: python

    import pandas as pd
    from snplib.format import make_ped

    input_data = pd.read_csv("file.txt")
    data_ped = make_ped(
        input_data, "SAMPLE_ID", "SNP", fid_col="SAMPLE_ID"
    )

    or

    data_ped = make_ped(
        input_data,
        "SAMPLE_ID",
        "SNP",
        fid_col="FAMILY_ID",
        father_col="father",
        mother_col="mother",
        sex_col="sex"
    )

Input data view::

   SAMPLE_ID          SNP
        1100  A A B B 0 0
        1101  A A B B B B
        1102  A A 0 0 B B
        1103  A A B B B B

    or

   SAMPLE_ID          SNP  FAMILY_ID  father  mother  sex
        1100  A A B B 0 0       1100       1       5    1
        1101  A A B B B B       1101       2       6    2
        1102  A A 0 0 B B       1102       3       7    1
        1103  A A B B B B       1103       4       8    0

Output data view::

    fid   sid father mother sex not_used          snp
   1100  1100      0      0   0        0  A A B B 0 0
   1101  1101      0      0   0        0  A A B B B B
   1102  1102      0      0   0        0  A A 0 0 B B
   1103  1103      0      0   0        0  A A B B B B

    or

    fid   sid father mother sex not_used          snp
   1100  1100      1      5   1        0  A A B B 0 0
   1101  1101      2      6   2        0  A A B B B B
   1102  1102      3      7   1        0  A A 0 0 B B
   1103  1103      4      8   0        0  A A B B B B


**fam** - https://www.cog-genomics.org/plink/1.9/formats#fam

.. code-block:: python

    import pandas as pd
    from snplib.format import make_fam

    input_data = pd.read_csv("file.txt", sep=" ")
    data_fam = make_fam(input_data, "SAMPLE_ID", "SAMPLE_ID")

    or

    make_fam(
        input_data,
        "SAMPLE_ID",
        "FAMILY_ID",
        father_col="father",
        mother_col="mother",
        sex_col="sex",
        pheno_col="pheno"
    )

Input data view::

   SAMPLE_ID  SNP
        1100  025
        1101  022
        1102  052
        1103  022

    or

   SAMPLE_ID  SNP  FAMILY_ID  father  mother  sex  pheno
       1100  025       1100       1       5    1     12
       1101  022       1101       2       6    2     13
       1102  052       1102       3       7    1     14
       1103  022       1103       4       8    0     15

Output data view::

     fid   sid father mother sex pheno
    1100  1100      0      0   0    -9
    1101  1101      0      0   0    -9
    1102  1102      0      0   0    -9
    1103  1103      0      0   0    -9

    or

     fid   sid father mother sex pheno
    1100  1100      1      5   1    12
    1101  1101      2      6   2    13
    1102  1102      3      7   1    14
    1103  1103      4      8   0    15


**lgen** - https://www.cog-genomics.org/plink/1.9/formats#lgen

.. code-block:: python

    import pandas as pd
    from snplib.format import make_lgen

    input_data = pd.read_csv("file.txt", sep=" ")
    data_lgen = make_lgen(
        input_data, "Sample ID", "SNP Name", ["Allele1 - AB", "Allele2 - AB"]
    )

Input data view::

     "SNP Name" "Sample ID" "Allele1 - AB" "Allele2 - AB" "GC Score" "GT Score"
                  ABCA12 107232207 A A 0.4048 0.8164
      ARS-BFGL-BAC-13031 107232207 B B 0.9083 0.8712
      ARS-BFGL-BAC-13039 107232207 A A 0.9005 0.9096
      ARS-BFGL-BAC-13049 107232207 A B 0.9295 0.8926
        ...
                   ABCA12 107237284 A A 0.4048 0.8164
       ARS-BFGL-BAC-13031 107237284 A B 0.9566 0.9257
       ARS-BFGL-BAC-13039 107237284 A B 0.3098 0.8555
       ARS-BFGL-BAC-13049 107237284 A B 0.8613 0.8319
        ...


Output data view::

    fid       sid            snp_name allele1 allele2
     1  107232207              ABCA12       A       A
     1  107232207  ARS-BFGL-BAC-13031       B       B
     1  107232207  ARS-BFGL-BAC-13039       A       A
     1  107232207  ARS-BFGL-BAC-13049       A       B
     1  107232207  ARS-BFGL-BAC-13059       A       B

     ...

     1  107237284              ABCA12       A       A
     1  107237284  ARS-BFGL-BAC-13031       A       B
     1  107237284  ARS-BFGL-BAC-13039       A       B
     1  107237284  ARS-BFGL-BAC-13049       A       B
     1  107237284  ARS-BFGL-BAC-13059       A       A
     ...



Statistics
----------

Poor quality or uninformative SNPs can be excluded from the analysis. This
helps to reduce noise and improve the accuracy of the results.


Call Rate
_________

The call rate for a given SNP is defined as the proportion of
individuals in the study for which the corresponding SNP information is
not missing. In the following example, we filter using a call rate of 95%,
meaning we retain SNPs for which there is less than 5% missing data.

**call rate marker**

Of the say, 54K markers in the chip, 50K have been genotyped for a
particular animal, the “call rate animal” is 50K/54K=93%

in_data::

        SNP_NAME SAMPLE_ID SNP
                ABCA12 1100 0
                 APAF1 1100 2
    ARS-BFGL-BAC-10172 1100 5
                ABCA12 1101 0
                 APAF1 1101 2
    ARS-BFGL-BAC-10172 1101 2
                ABCA12 1102 0
                 APAF1 1102 5
    ARS-BFGL-BAC-10172 1102 2
                ABCA12 1103 0
                 APAF1 1103 2
    ARS-BFGL-BAC-10172 1103 2
                ABCA12 1104 5
                 APAF1 1104 1
    ARS-BFGL-BAC-10172 1104 1
                ABCA12 1105 0
                 APAF1 1105 2
    ARS-BFGL-BAC-10172 1105 5
                ABCA12 1106 0
                 APAF1 1106 1
    ARS-BFGL-BAC-10172 1106 2
                ABCA12 1107 5
                 APAF1 1107 2
    ARS-BFGL-BAC-10172 1107 1
                ABCA12 1108 0
                 APAF1 1108 2
    ARS-BFGL-BAC-10172 1108 2
                ABCA12 1109 0
                 APAF1 1109 2
    ARS-BFGL-BAC-10172 1109 2
                ABCA12 1110 5
                 APAF1 1110 2
    ARS-BFGL-BAC-10172 1110 2

.. code-block:: python

    import pandas as pd
    from snplib.statistics import call_rate

    input_data = pd.read_csv("file.txt", sep=" ")
    result = call_rate(data=input_data, id_col="SNP_NAME", snp_col="SNP")

result::

                 SNP_NAME       SNP
                   ABCA12  0.727273
                    APAF1  0.909091
       ARS-BFGL-BAC-10172  0.818182

**call rate animal**

Of the say, 900 animals genotyped for marker CL635944_160.1, how many
have actually been successfully read? Assume that 600 have been read, then
the “call rate marker” is 600/900 = 67%

in_data::

                  SNP_NAME SAMPLE_ID SNP
                    ABCA12     14814   0
        ARS-BFGL-BAC-13031     14814   2
        ARS-BFGL-BAC-13039     14814   0
        ARS-BFGL-BAC-13049     14814   1
        ARS-BFGL-BAC-13059     14814   1
        ARS-BFGL-BAC-13086     14814   0
        ARS-BFGL-BAC-13093     14814   1
        ARS-BFGL-BAC-13110     14814   5
        ARS-BFGL-BAC-13111     14814   0
        ARS-BFGL-BAC-13113     14814   1
        ARS-BFGL-BAC-15633     14814   0
        ARS-BFGL-BAC-15634     14814   0
        ARS-BFGL-BAC-15637     14814   0
        ARS-BFGL-BAC-15659     14814   0
        ARS-BFGL-BAC-15668     14814   5
        ARS-BFGL-BAC-15708     14814   0
        ARS-BFGL-BAC-15718     14814   0
                    ABCA12     14815   0
        ARS-BFGL-BAC-13031     14815   1
        ARS-BFGL-BAC-13039     14815   1
        ARS-BFGL-BAC-13049     14815   1
        ARS-BFGL-BAC-13059     14815   0
        ARS-BFGL-BAC-13086     14815   1
        ARS-BFGL-BAC-13093     14815   5
        ARS-BFGL-BAC-13110     14815   2
        ARS-BFGL-BAC-13111     14815   1
        ARS-BFGL-BAC-13113     14815   2
        ARS-BFGL-BAC-15633     14815   0
        ARS-BFGL-BAC-15634     14815   2
        ARS-BFGL-BAC-15637     14815   2
        ARS-BFGL-BAC-15659     14815   2
        ARS-BFGL-BAC-15668     14815   5
        ARS-BFGL-BAC-15708     14815   1
        ARS-BFGL-BAC-15718     14815   2

.. code-block:: python

    import pandas as pd
    from snplib.statistics import call_rate

    input_data = pd.read_csv("file.txt", sep=" ")
    result = call_rate(data=data_df, id_col="SAMPLE_ID", snp_col="SNP")

result::

      SAMPLE_ID       SNP
          14814  0.882353
          14815  0.882353


Frequence Allele
________________

The allele frequency represents the incidence of a gene variant in a
population.

**allele freq**

in_data::

        SNP_NAME SAMPLE_ID SNP
                ABCA12 1100 0
                 APAF1 1100 2
    ARS-BFGL-BAC-10172 1100 5
                ABCA12 1101 0
                 APAF1 1101 2
    ARS-BFGL-BAC-10172 1101 2
                ABCA12 1102 0
                 APAF1 1102 5
    ARS-BFGL-BAC-10172 1102 2
                ABCA12 1103 0
                 APAF1 1103 2
    ARS-BFGL-BAC-10172 1103 2
                ABCA12 1104 5
                 APAF1 1104 1
    ARS-BFGL-BAC-10172 1104 1
                ABCA12 1105 0
                 APAF1 1105 2
    ARS-BFGL-BAC-10172 1105 5
                ABCA12 1106 0
                 APAF1 1106 1
    ARS-BFGL-BAC-10172 1106 2
                ABCA12 1107 5
                 APAF1 1107 2
    ARS-BFGL-BAC-10172 1107 1
                ABCA12 1108 0
                 APAF1 1108 2
    ARS-BFGL-BAC-10172 1108 2
                ABCA12 1109 0
                 APAF1 1109 2
    ARS-BFGL-BAC-10172 1109 2
                ABCA12 1110 5
                 APAF1 1110 2
    ARS-BFGL-BAC-10172 1110 2

.. code-block:: python

    import pandas as pd
    from snplib.statistics import allele_freq

    input_data = pd.read_csv("file.txt", sep=" ")
    result = allele_freq(data=input_data, id_col="SNP_NAME", seq_col="SNP")

result::

                 SNP_NAME    SNP
                   ABCA12  0.000
                    APAF1  0.900
       ARS-BFGL-BAC-10172  0.889

The minor allele frequency is therefore the frequency at which the
minor allele occurs within a population.

**maf**

.. code-block:: python

    from snplib.statistics import minor_allele_freq as maf

    result = maf(0.22)  # result == 0.22


HWE (Hardy-Weinberg equilibrium)
________________________________

The Hardy-Weinberg equilibrium is a principle stating that the genetic
variation in a population will remain constant from one generation to the
next in the absence of disturbing factors.
https://www.nature.com/scitable/definition/hardy-weinberg-equilibrium-122/

To test the hypothesis that the data are within the HWE, a statistic a chi2
distribution with 1 degree of freedom:

.. code-block:: python

    from snplib.statistics import hwe_test

    result = hwe_test(seq_snp=pd.Series(list(map(int, "2212120"))), freq=0.714)  # True
    result = hwe_test(seq_snp=pd.Series(list(map(int, "02011015010000500"))), freq=0.2)  # True
    result = hwe_test(seq_snp=pd.Series(list(map(int, "000000000102"))), freq=0.125)  # False


The p-value used here is:

.. code-block:: python

    from snplib.statistics import hwe

    hom1 = 10
    hets = 500
    hom2 = 5000

    result = hwe(hets, hom1, hom2)  # 0.6515718999145375 (p-value)

Once the data have been prepared, statistical analysis to identify associations,
patterns, or relationships between SNPs and the phenotypes or diseases of
interest (GWAS). phenotypes or diseases of interest (GWAS).

Parentage
---------
https://www.icar.org/Documents/GenoEx/ICAR%20Guidelines%20for%20Parentage%20Verification%20and%20Parentage%20Discovery%20based%20on%20SNP.pdf


.. note::
    A list of isag verification and discovery macerators can be found here.
    See Appendix list - https://www.icar.org/Guidelines/04-DNA-Technology.pdf

    List of SNP to be used for either parentage verification or parentage discovery (appendix 11):
    https://www.icar.org/Guidelines/04-DNA-Technology-App-11-SNP-list-for-parentage-verification-or-discovery.pdf


Verification
____________

Verification of paternity according to ICAR recommendations.

input data::

                       SNP_Name      ID41988163  ID10512586
    0                    ABCA12               0           0
    1                     APAF1               2           2
    2        ARS-BFGL-BAC-10172               2           2
    3         ARS-BFGL-BAC-1020               1           1
    4        ARS-BFGL-BAC-10245               1           1
    ..                      ...             ...         ...
    239  Hapmap55441-rs29010990               1           1
    240  Hapmap59876-rs29018046               1           0
    241  Hapmap60017-rs29023471               2           1
    242           UA-IFASA-5034               0           1
    243           UA-IFASA-6532               0           0


.. code-block:: python

    from snplib.parentage import Verification, isag_verif

    input_data = pd.read_csv("file.txt", sep=" ")

    obj_verification = Verification(isag_marks=isag_verif().markers)
    result = obj_verification.check_on(
        data=input_data,
        descendant="ID41988163",
        parent="ID10512586",
        snp_name_col="SNP_Name"
    )

    # Result
    print(obj_verification.num_conflicts)  # 31
    print(obj_verification.status)  # "Excluded"


Discovery
_________

Search for paternity according to ICAR recommendations.

input data::

                   SNP_Name      ID41988163  ID10512586
    0                ABCA12               0           0
    1                 APAF1               2           2
    2    ARS-BFGL-BAC-10172               2           2
    3     ARS-BFGL-BAC-1020               1           1
    4    ARS-BFGL-BAC-10245               1           1
    ..                  ...             ...         ...
    617       UA-IFASA-5034               0           1
    618       UA-IFASA-6154               2           0
    619       UA-IFASA-6532               0           0
    620       UA-IFASA-8658               1           0
    621       UA-IFASA-8833               0           0

.. code-block:: python

    from snplib.parentage import Discovery, isag_disc

    input_data = pd.read_csv("file.txt", sep=" ")

    obj_discovery = Discovery(isag_marks=isag_disc().markers)
    result = obj_discovery.search_parent(
        data=input_data,
        descendant="ID41988163",
        parents="ID10512586",
        snp_name_col="SNP_Name"
    )

    # Result
    print(obj_discovery.num_conflicts)  # 77
    print(obj_discovery.status)  # "Excluded"
    print(obj_discovery.perc_conflicts)  # 14.86 %
