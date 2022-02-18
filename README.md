# MetaLP
MetaLP is a protein inference algorithm for shotgun proteomics analysis of microbial communities. Two key innovations in MetaLP were the integration of taxonomic abundances as prior information and the formulation of protein inference as a linear programming problem. It was optimized for metaproteomics to address the degenerate peptides and ‘one-hit wonders’ issues. The objective of MetaLP is to produce substantially higher numbers of protein identifications in complex metaproteomics datasets.

## Setup and installation
### Dependency
* Anaconda3 for python environment
* numpy >= 1.19.5
* gurobipy == 9.1.2
* Note that python package for gurobi requires a license. Gurobi provides free academic licenses after the registration (https://pages.gurobi.com/registration). Detail instruction for the license registration is available (https://www.gurobi.com/downloads/end-user-license-agreement-academic/).

## Quick Start for MetaLP
```
python LP2_species_V2.py -i input/P2.identification_1%.txt -o P2_result.txt -g P2_filteringResult.txt -d pro2otu/P2_ProToOTU.tsv -f 0.01 -p probability/protein_probability_P2.txt

```
* MetaLP requires three input files shown as below:
```
-i: a file contains the peptide identification scores for each identified peptide and corresponding parent proteins
-d: a file contains the relationship between protein identifier with the OTUs which the protein is belonging to
-p: a file contains OTUs (operational Taxonomic Units) and their estimated probability
```
* MetaLP generates two output files shown as below:
```
-o: protein inference result, it contains protein groups and their score inferred by MetaLP
-g: filtering result by assigned FDR (False Discovery Rate). The default FDR is 0.01, which can be passed by the option -f
```
* The detail user manual will show if using:
```
python LP2_species_V2.py -h
```

## Metagenomics analysis commands
### Dependency
* MetaWrap pipeline (https://github.com/bxlab/metaWRAP)

### Assembling and binning
```
conda activate metaWRAP-env
metawrap assembly -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -m 200 -t 30 --metaspades -o ASSEMBLY

bowtie2-build -f /work/fs0199/LP_human_gut/ASSEMBLY/final_genomics_assembly.fasta final --threads 32
bowtie2 -1 SRR3313113_2_val_2.fq -2 SRR3313113_2_val_2.fq -p 32 -x final -S final.sam

samtools view -@ 16 -b -S final.sam -o final.bam
samtools sort -@ 32 -l 9 -O BAM final.bam -o final.sorted.bam

jgi_summarize_bam_contig_depths --outputDepth final.depth.txt final.sorted.bam
metabat2 -t 32 -i /work/fs0199/LP_human_gut/ASSEMBLY/final_genomics_assembly.fasta -a final.depth.txt -o output.human_gut
```
## Utils
### Utils for OTU probability estimation
* In the folder OTU_probability_utils we provide our tools to calculate the OTU probability. The file MRC.out is to calculate OTUs' probability for soil and marine metaproteome (unknown OTU) and the file MockMRC.out is to calculate OTUs' probability for the mock community (known OTU). A sample command is shown as below:
```
./MRC.out -i samfile.sam -c OTU_reference_genomes_directory
./MockMRC.out -i samfile.sam -o outputFile
```
* samfile is the sequence alignment mapping file. It is output file by BBsplit for mock community and by Bowtie2 for soil/marine microbial communities seperately. OTU reference geonomes directory is the binning results by metaBAT2
* human_gut_16s_count.txt is the reference file to calculated OTU probability for human gut. It is organized by result from vsearch and Blast

### Utils for PeptideProphet
* The script to parse PeptideProphet result from  .xml format to tab-delimited format is provided as "peptideprophetParser.py". The user manual is available by:
```
python peptideprophetParser.py -h
```
