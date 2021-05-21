# TSEBRA: Transcript Selector for BRAKER
### Introduction
TSEBRA is a combiner tool that selects transcripts from gene predictions based on the support by extrisic evidence in form of introns and start/stop codons. It was developed to combine BRAKER1<sup name="a1">[1](#ref1)</sup> and BRAKER2<sup name="a2">[2](#ref2)</sup> predicitons to increase their accuracies.

## Prerequisites
Python 3.5.2 or higher is required.

## Installation
Download TSEBRA:
```console
git clone https://github.com/Gaius-Augustus/TSEBRA
```
Or download TSEBRA as submodule of BRAKER with:
```console
git clone --recurse-submodules https://github.com/Gaius-Augustus/BRAKER
```

## Usage
The main script is ```./bin/tsebra.py```. For usage information run ```./bin/tsebra.py --help```.

## Input Files
TSEBRA takes a list of gene prediciton files, a list of hintfiles and a configuration file as mandatory input.

#### Gene Predictions 
The gene prediction files have to be in gtf format. This is the standard output format of a BRAKER or AUGUSTUS<sup name="a3">[3,](#ref3)</sup><sup name="a4">[4](#ref4)</sup> gene prediciton.

Example:
```console
2L      AUGUSTUS        gene    83268   87026   0.88    -       .       g5332
2L      AUGUSTUS        transcript      83268   87026   0.88    -       .       g5332.t1
2L      AUGUSTUS        intron  84278   87019   1       -       .       transcript_id "file_1_file_1_g5332.t1"; gene_id "file_1_file_1_g5332";
2L      AUGUSTUS        CDS     87020   87026   0.88    -       0       transcript_id "file_1_file_1_g5332.t1"; gene_id "file_1_file_1_g5332";
2L      AUGUSTUS        exon    87020   87026   .       -       .       transcript_id "file_1_file_1_g5332.t1"; gene_id "file_1_file_1_g5332";
```

#### Hint Files
The hints files have to be in gff format, the last column must include an attribute for the source for the hint with 'src=' and can include the number of hints supporting the gene structure segment with 'mult='. This is the standard file format of the ```hintfiles.gff``` in a BRAKER working directory.

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
1. The weight for each hint source. (A weight is set to 1, if the weight for a source is not determined in the cfg file.)
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
The recommended and most common usage for TSEBRA is to combine the resulting ```braker.gtf``` files of a BRAKER1 and a BRAKER2 run using the hintsfile.gff from both working directories. However, TSEBRA can be applied to any number (>1) of gene predictions and hint files as long as they are in the correct format. 

A common case might be that a user wants to annotate a novel genome with BRAKER and has:
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
2. Make sure that the gene and transcript IDs of the gene prediction files are in order (this step is optional)
```console
./bin/fix_gtf_ids.py --gtf braker1_out/braker.gtf --out braker1_fixed.gtf
./bin/fix_gtf_ids.py --gtf braker2_out/braker.gtf --out braker2_fixed.gtf
```
3. Combine predicitons with TSEBRA
```console
./bin/tsebra.py -g braker1_fixed.gtf,braker2_fixed.gtf -c default.cfg \ 
    -e braker1_out/hintsfile.gff,braker2_out/hintsfile.gff \
    -o braker1+2_combined.gtf
```
The combined gene prediciton is ```braker1+2_combined.gtf```.

## Example
A small example is located at ```example/```. Run ```./example/run_prevco_example.sh``` to execute the example and to check if TSEBRA runs properly. 

## Licence
All source code, i.e. `bin/*.py` are under the [Artistic License](bin/LICENSE.txt) (see <https://opensource.org/licenses/Artistic-2.0>).

## References
<b id="ref1">[1]</b> Hoff, Katharina J, Simone Lange, Alexandre Lomsadze, Mark Borodovsky, and Mario Stanke. 2015. “BRAKER1: Unsupervised Rna-Seq-Based Genome Annotation with Genemark-et and Augustus.” *Bioinformatics* 32 (5). Oxford University Press: 767--69.[↑](#a1)

<b id="ref2">[2]</b> Tomas Bruna, Katharina J. Hoff, Alexandre Lomsadze, Mario Stanke and Mark Borodvsky. 2021. “BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database." *NAR Genomics and Bioinformatics* 3(1):lqaa108.[↑](#a2)

<b id="ref3">[3]</b> Stanke, Mario, Mark Diekhans, Robert Baertsch, and David Haussler. 2008. “Using Native and Syntenically Mapped cDNA Alignments to Improve de Novo Gene Finding.” *Bioinformatics* 24 (5). Oxford University Press: 637--44.[↑](#a3)

<b id="ref4">[4]</b> Stanke, Mario, Oliver Schöffmann, Burkhard Morgenstern, and Stephan Waack. 2006. “Gene Prediction in Eukaryotes with a Generalized Hidden Markov Model That Uses Hints from External Sources.” *BMC Bioinformatics* 7 (1). BioMed Central: 62.[↑](#a4)
