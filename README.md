# <ins>PrEvCo</ins>: Gene <ins>Pr</ins>edcition and Extrinsic <ins>Ev</ins>idence <ins>Co</ins>mbiner Tool
### Introduction
PrEvCo is a combiner tool that combines gene predictions based on the support by extrisic evidence in form of introns and start/stop codons. It was developed to combine BRAKER1 and BRAKER2 predicitons to increase their accuracy.

## Prerequisites
Python 3.5.2 or higher is required.

## Installation
Download PrEvCo with 
```console
git clone ToDo: ENTER LINK 
```

## Usage
The main script is ```./bin/prevco.py```. For usage information run ```./bin/prevco.py --help```.

## Input Files
PrEvCo needs a list of gene prediciton files, a list of hintfiles and a configuration file as input.

#### Gene Predictions 
The gene prediction files needs to be in gtf format. This is the standard output format of a BRAKER or AUGUSTUS gene prediciton.

Example:
```console
2L      AUGUSTUS        gene    83268   87026   0.88    -       .       g5332
2L      AUGUSTUS        transcript      83268   87026   0.88    -       .       g5332.t1
2L      AUGUSTUS        intron  84278   87019   1       -       .       transcript_id "file_1_file_1_g5332.t1"; gene_id "file_1_file_1_g5332";
2L      AUGUSTUS        CDS     87020   87026   0.88    -       0       transcript_id "file_1_file_1_g5332.t1"; gene_id "file_1_file_1_g5332";
2L      AUGUSTUS        exon    87020   87026   .       -       .       transcript_id "file_1_file_1_g5332.t1"; gene_id "file_1_file_1_g5332";
```

#### Hint Files
The hint files have to be in gff format, the last column must include an attribute for the source for the hint with 'src=' and can include the number of hints supporting the gene structure segment with 'mult='. This is the standard file format of the ```hintfiles.gff``` in a BRAKER working directory.

Example:
```console
2L      ProtHint        intron  279806  279869  2       +       .       src=P;mult=25;pri=4;al_score=0.437399;
2L      ProtHint        intron  275252  275318  2       -       .       src=P;mult=19;pri=4;al_score=0.430006;
2L      ProtHint        stop    293000  293002  1       +       0       grp=7220_0:002b08_g42;src=C;pri=4;
2L      ProtHint        intron  207632  207710  1       +       .       grp=7220_0:002afa_g26;src=C;pri=4;
2L      ProtHint        start   207512  207514  1       +       0       grp=7220_0:002afa_g26;src=C;pri=4;
```

#### Configuration File
The configuration file has to include three types of parameter:
1. The weight for each hint source. A weight is set to 1, if the weight for a source is not determined in the cfg file.
2. Required fraction of supported introns or supported start/stop-codons for a transcript.
3. Allowed difference between two overlapping transcripts for each feature type.

Example:
```console
# src weights
P 0.1
E 10
C 5
M 1
# Low evidence support
intron_support 0.75
stasto_support 1
# Feature differences 
e_1 0
e_2 0.5
e_3 25
e_4 10
```


## Use Case
The recommended and most common usage for PrEvCo is to combine the resultingbraker.gtffiles of a BRAKER1 and a BRAKER2 run using thehintsfile.gff from both working directories. However, PrEvCo can be applied to any number (>1) of gene predictions and hint files as long as they are in the correct format. A common case might be that a user wants to annotate a novel genome with BRAKER and has:
A typical case for running BRAKER and PrEvCo would be, if you have
* a novel genome with repeats masked: ```genome.fasta.masked```,
* hints for intron positions from RNA-seq reads```rna_seq_hints.gff```,
* database of homologous proteins: ```proteins.fa```.

1. Run BRAKER1 and BRAKER2 for example with
```console
### BRAKER1
braker.pl --genome=genome.fasta.masked --hints=rna_seq_hints.gff \ 
            --softmasking --species=species_name --workingdir=braker1_out
    
### BRAKER2
braker.pl --genome=genome.fasta.masked --prot_seq=proteins.fa \ 
    --softmasking --species=species_name --epmode --prg=ph \ 
    --workingdir=braker2_out
```
2. Make sure that the gene and transcript IDs of the gene prediction files are in order (This step is optional)
```console
./bin/fix_gtf_ids.py --gtf braker1_out/braker.gtf --out braker1_fixed.gtf
./bin/fix_gtf_ids.py --gtf braker2_out/braker.gtf --out braker2_fixed.gtf
```
3. Combine predicitons with PrEvCo
```console
./bin/prevco.py -g braker1_fixed.gtf,braker2_fixed.gtf -c default.cfg \ 
    -e braker1_out/hintsfile.gff,braker2_out/hintsfile.gff \
    -o braker1+2_combined.gtf
```
The combined gene prediciton is ```braker1+2_combined.gtf```.

## Example
A small example is located at ```example/```. Run ```./example/run_prevco_example.sh``` to execute the example and to check if PrEvCo runs properly. 

## Reference
ToDo

