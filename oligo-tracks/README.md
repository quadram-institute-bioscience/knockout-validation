# Map oligos to reference

## Dependencies

The required tools can be installed via Miniconda:

```
conda install -y bwa samtools covtobed
```

## Usage

The program requires a reference and a set of oligos, both in FASTA format:

```
perl align_oligos.pl -r {REFERENCE.fa} -i {OLIGOS.fa} -o {OUTPUT_NAME}
```

The program will create:
* {OUTPUT_NAME}.bam (aligned oligos, intermediate file)
* {OUTPUT_NAME}.bed (track with mapped oligos)