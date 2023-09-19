# plate-demux

![example workflow](https://github.com/gartician/plate-demux/actions/workflows/test.yaml/badge.svg)
[![python minimum](https://img.shields.io/badge/python->=3-<COLOR>.svg)](https://shields.io/)
![Maintainer](https://img.shields.io/badge/maintainer-gartician-blue)

Demultiplex combinatorial indexing sequencing data (e.g. S3-ATAC-Seq) performed in 96-well plates.

# Usage

```
usage: plate-demux.py [-h] -R1 R1 -R2 R2 -o OUTDIR -c CONFIG [-b BUFFER] [-v]

plate-demux.py demultiplexes ATAC-Seq/S3-ATAC-Seq reads performed in 96-well plates.
This script assumes that a Tn5 index (e.g. an 8bp DNA barcode) tags for a specific
biological sample, and that index appears in a specific position in a 96-well plate
(e.g. index ACTAAGTAA in A12). By relating the index to the coordinates of a plate,
this sript will demultiplex a FASTQ file of mixed samples into separate files.

This script accepts paired FASTQ files processed with unidex () and a configuration file.
The FASTQ file should contain a mix of samples, while the configuration is a table
that specifies an index and its coordinates in a 96-well plate. The output is a folder
that contains foward and reverse reads per individual sample. Due to potentially large
input FASTQ files, plate-demux.py allows users to process the data piece-wise by
loading <n> reads into memory at a time via the --buffer argument.

optional arguments:
  -h, --help            show this help message and exit
  -R1 R1                Forward read (R1) of a paired sequencing run.
  -R2 R2                Reverse read (R2) of a paired sequencing run.
  -o OUTDIR, --outdir OUTDIR
                        Output directory with forward and reverse reads per sample.
  -c CONFIG, --config CONFIG
                        File containing tn5 indices, well positions, and biological samples.
  -b BUFFER, --buffer BUFFER
                        Maximum number of reads to load into memory. The program generally
                        requires approximately 1GB of memory per 1M reads.
  -v, --verbose         Adjust the verbosity of the program.
```

# Quickstart

```bash
# bulk demux = load input FASTQs into memory
python plate-demux.py \
  -R1 test-data/mixed-samples-R1.fastq.gz \
  -R2 test-data/mixed-samples-R2.fastq.gz \
  -c config.txt \
  -o output-data \
  --verbose

# piece-wise demux = load <n> reads into memory, export, and repeat.
python plate-demux.py \
  -R1 test-data/mixed-samples-R1.fastq.gz \
  -R2 test-data/mixed-samples-R2.fastq.gz \
  -c config.txt \
  -o output-data \
  -b 100000 \
  --verbose
```

# Input

## FASTQ file

Input FASTQ files should be the output of [unidex](https://github.com/ohsu-cedar-comp-hub/unidex) and contain a mix of samples. In the example files in the `test-data` folder, these FASTQ files should have the Tn5 barcode as the last 8 base pairs of a FASTQ record's header such as the example below:

```
@TAAGGCGATTAACTCATAAAGTCCAA:1#0/1
GTGATTGGTGGGTCATTATGTGTTGTCGTGCAGGTAGAGGCCTGTCTCTTATACACATCTATTGGACTTAGATCGGAAGAGCACACGTCT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@AGGCAGAACTATGGCCTAACCTTGGC:2#0/1
TATTCATCCTATGTGGGTAATTGAGGAGTATGCTAAGATTTTGCGTAGCTGGGTTTGGTTTAATCCACCTCAACTGCCTGCTATCTGTCT
+
CCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCC-C--CC;CCCC;CC;CCCCCCCCCCCCCCCCCCC;CCCC-;CCCCCCCC;CCCCCCC
@TCCTGAGCCTATGGCCTACCTTCACC:3#0/1
GCTAGAGACTAAGGGAGGGGCAGTTTCAGAGAACCGAGGGGGGAAAAGTCTGTGTTTCTTTAGTCTCACATAGTGTCTTGAATAACCAGA
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCC;CCCCCCCCCCCCCCCCCCC
```

Where AAGTCCAA, ACCTTGGC, and CCTTCACC (from records 1, 2, and 3) correspond to a position in a 96-well plate specified in the configuration file.

NOTE: Due to potentially large input FASTQ files, plate-demux.py allows users to process the data piece-wise by loading `--buffer <n>` reads into memory at a time via the --buffer argument. This program generally requires about 1GB of memory per 1M reads.

## Configuration file

The configuration file `config.txt` must contain a "Tn5_index" and "sample" column. The program will associate an index to a sample, and a sample can use multiple Tn5 indices. In the example files, two biological samples were sequenced where wildtype-1 utilized columns 1-6 and mutant-1 were in columns 7-12. While the well positions are not explicitly used to demultiplex samples, the `config.txt` file provides a clear annotation for bench and computational scientists to unambiguously identify samples.

```
Well   Row     Column  Name    Truseq_R2_SBS12_partial Tn5_index       sample
D1      D       1       SBS12_18_UME_sci_37     CGTGTGCTCTTCCGATCT      ATGTAAGT        wildtype-1
D2      D       2       SBS12_18_UME_sci_38     CGTGTGCTCTTCCGATCT      GCACGGAC        wildtype-1
D3      D       3       SBS12_18_UME_sci_39     CGTGTGCTCTTCCGATCT      GGTACCTT        wildtype-1
D4      D       4       SBS12_18_UME_sci_40     CGTGTGCTCTTCCGATCT      AACGTTCC        wildtype-1
D5      D       5       SBS12_18_UME_sci_41     CGTGTGCTCTTCCGATCT      GCAGAATT        wildtype-1
D6      D       6       SBS12_18_UME_sci_42     CGTGTGCTCTTCCGATCT      ATGAGGCC        wildtype-1
D7      D       7       SBS12_18_UME_sci_43     CGTGTGCTCTTCCGATCT      ACTAAGAT        mutant-1
D8      D       8       SBS12_18_UME_sci_44     CGTGTGCTCTTCCGATCT      GTCGGAGC        mutant-1
D9      D       9       SBS12_18_UME_sci_45     CGTGTGCTCTTCCGATCT      CCGCGGTT        mutant-1
D10     D       10      SBS12_18_UME_sci_46     CGTGTGCTCTTCCGATCT      TTATAACC        mutant-1
D11     D       11      SBS12_18_UME_sci_47     CGTGTGCTCTTCCGATCT      GGACTTGG        mutant-1
D12     D       12      SBS12_18_UME_sci_48     CGTGTGCTCTTCCGATCT      AAGTCCAA        mutant-1
```

# Output

The output of this program is a folder that contains forward and reverse reads for every biological sample specified in the `config.txt`. Below contains the example output of this script:

```bash
$ python plate-demux.py \
  -R1 test-data/mixed-samples-R1.fastq.gz \
  -R2 test-data/mixed-samples-R2.fastq.gz \
  -c config.txt \
  -o output-data

$ tree output-data
output-data/
├── mutant-1_R1.fastq.gz
├── mutant-1_R2.fastq.gz
├── wildtype-1_R1.fastq.gz
└── wildtype-1_R2.fastq.gz
```

# References

Mulqueen, R.M., Pokholok, D., O’Connell, B.L. et al. High-content single-cell combinatorial indexing. Nat Biotechnol 39, 1574–1580 (2021). https://doi.org/10.1038/s41587-021-00962-z

Vitak, S., Torkenczy, K., Rosenkrantz, J. et al. Sequencing thousands of single-cell genomes with combinatorial indexing. Nat Methods 14, 302–308 (2017). https://doi.org/10.1038/nmeth.4154

[Advances in single-cell combinatorial indexing](https://www.takarabio.com/learning-centers/automation-systems/icell8-introduction/advances-in-single-cell-combinatorial-indexing)
