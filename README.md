# PrEvCo: gene Predcition and extrinsic Evidence Combiner
### Introduction
PrEvCo is a combiner tool that combines gene predictions based on the support by extrisic evidence in form of introns and start/stop codons. It was developed to combine BRAKER1 and BRAKER2 predicitons to increase their accuracy.

### Prerequisites
Python 3.5.2 or higher is required.

### Installation
Download PrEvCo with 
```console
git clone ToDo: ENTER LINK 
```

### Usage
The main script is ```./bin/prevco.py```. For usage information run ```./bin/prevco.py --help```.

### Input Files
#### Gene Predictions 
ToDo
#### Hint Files
ToDo
#### Configuration File
ToDo
### Use Case
PrEvCo is intended to be used together with BRAKER for creating a BRAKER-only prediciton supported by RNA-seq and homologous protein data.
A typical case for running BRAKER and PrEvCo would be, if you have
* ```genome.fasta.masked```: a novel genome masked for repeats,
* ```rna_seq_hints.gff```: intron position evidence from RNA-seq reads,
* ```proteins.tab```: database of homologous proteins.
1. Run BRAKER1 and BRAKER2 for example with
```console
### BRAKER1
braker.pl --genome=genome.fasta.masked --hints=rna_seq_hints.gff \ 
            --softmasking --species=species_name --workingdir=braker1_out
    
### BRAKER2
braker.pl --genome=genome.fasta.masked --prot_seq=proteins.tab \ 
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

### Example
A small example is located at ```example/```. Run ```./example/run_prevco_example.sh``` to execute the example check if PrEvCo runs properly. 

### Reference
ToDo

