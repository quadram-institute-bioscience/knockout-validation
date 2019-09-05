## Scripts

 - *check_knocked_out.py*: script to parse a count matrix generated with featureCounts (subread package), reporting the genes that could be KO in a sample.
 - *grep_fasta.pl*: detect a motif in a multifasta sequence and report position and flanking sequences
 
### Output

The list of knocked out genes has been produced with the command:
```
featureCounts -a ../ref/CP009273_annotation.gtf -o ../coverage/counts.txt -t CDS aln/*.bam
python detect.py -i ../coverage/counts.txt -a ../ref/CP009273_annotation.gff  -n 100 
```
