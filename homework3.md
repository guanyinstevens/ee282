# Commands Used in Analysis and Corresponding Outputs
## Downloading data
`wget "http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz" -O dmel-all-chromosome-r6.48.fasta.gz`
This downloads the data from the internet, reassigning it to the file named dmel-all-chromosome-r6.48.fasta.gz
## Checking the quality of the data performing a checksum, using the md5 checksum since the data is originally from FlyBase
`md5sum dmel-all-chromosome-r6.48.fasta.gz`
Can open the md5sum.txt available and check that the output matches the line for the specific file. 
`wget "http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/md5sum.txt" -O md5sum.txt`
`md5sum -c md5sum.txt`
Outputs what file you have downloaded, and if the data is OK. 
An additional command: `md5sum *.gz > md5check.txt` allows us to check again, while creating our own md5 check file. 
## Processing and analyzing the FASTA file
`faSize dmel-all-chromosome-r6.48.fasta.gz`
This summarizes the FASTA file, listing how many total bases, N's, as well as the size and standard deviation of the file. 
`faSize -detailed dmel-all-chromosome-r6.48.fasta.gz`
The output will show the chromosomes, mapped scaffolds, unmapped scaffold, and rDNA. 
## Counting the number of nucleotides
`zcat dmel-all-chromosome-r6.48.fasta.gz | grep -v '^>' | tr -d '\n' | wc -c > nucleotides.txt`
`less nucleotides`
Will open file nucleotides with the number 143726002. This allows to remove the header/geneID by using grep. The command `tr -d '\n'` will delete any spacing in the line, essentially turning into one line. We can then count the characters of that line using wc -c, telling us the number of nucleotides. 
## Counting the number of sequences
`zcat dmel-all-chromosome-r6.48.fasta.gz | grep -c '^>' > sequences.txt`
`less sequences`
Using grep to count the number of things that starts with ^> allows us to count how many gene sequences are present. 
Will open the file sequences that will tell us that we have 1870 sequences
## Counting the number of N's
`zcat dmel-all-chromosome-r6.48.fasta.gz | grep -v '^>' | tr -d '\n' | grep -o N | wc -l > N.txt`
`less N.txt`
1152978. This summarizes that there are 1152978 unknown nucleotides.
## downloading the annotated gtf file
`wget "http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/dmel-all-r6.48.gtf.gz" -O dmel-all-r6.48.gtf.gz`
## checking data quality of the gtf file
`md5sum dmel-all-r6.48.gtf.gz`
Can open the md5sum.txt available and check that the output matches the line for the specific file. 
`md5sum *gtf.gz > md5check.gtf.txt`
`less md5check.gtf.txt`
Will output file that we have downloaded, containing the checksum. Checking that they are the same will tell us that our data is ok. 
## Sorting features from most to least common
`bioawk -c gff '{print $3}' dmel-all-r6.48.gtf.gz | sort | uniq -c | sort -nr > sorted.features.txt`
`less sorted.features.txt`
The output should be:
190050 exon
 163242 CDS
  46802 5UTR
  33738 3UTR
  30885 start_codon
  30825 stop_codon
  30799 mRNA
  17896 gene
   3053 ncRNA
    485 miRNA
    365 pseudogene
    312 tRNA
    300 snoRNA
    262 pre_miRNA
    115 rRNA
     32 snRNA
## Total number of genes on the X chromosome (Individually)
`bioawk -c gff '{if ( $1 == "X" && $3 == "gene" ) print $0}' dmel-all-r6.48.gtf.gz | sort | uniq -c | wc -l > X.genes.txt`
`less X.genes.txt`
2708. There are 2708 genes on the X chromosome.
## Total number of genes on the 2L arm (Individually)
`bioawk -c gff '{if ( $1 == "2L" && $3 == "gene" ) print $0}' dmel-all-r6.48.gtf.gz | sort | uniq -c | wc -l > 2L.genes.txt`
`less 2L.genes.txt`
3515
While we could continue to do that individually for each arm, we can sort them all at once. 
## Summarizing the number of genes per chromosome arm.
We can do this using a similar method to when we sorted the number of features. 
`bioawk -c gff '$3 == "gene" {print $1}' dmel-all-r6.48.gtf.gz | sort | uniq -c | sort -nr > Sorted.Chromosomes.txt`

`less Sorted.Chromosomes.txt`
4227 3R
   3653 2R
   3515 2L
   3489 3L
   2708 X
    114 4
    113 Y
     38 mitochondrion_genome
     21 rDNA
      2 Unmapped_Scaffold_8_D1580_D1567
      2 211000022280494
      1 211000022280703
      1 211000022280481
      1 211000022280347
      1 211000022280341
      1 211000022280328
      1 211000022279681
      1 211000022279392
      1 211000022279264
      1 211000022279188
      1 211000022279165
      1 211000022278760
      1 211000022278449
      1 211000022278436
      1 211000022278279
