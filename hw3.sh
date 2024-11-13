## Commands for Homework 3
# Summaries of Genome
`wget "http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz" -O dmel-all-chromosome-r6.48.fasta.gz`
`md5sum dmel-all-chromosome-r6.48.fasta.gz`
`wget "http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/md5sum.txt" -O md5sum.txt`
`md5sum -c md5sum.txt`
`md5sum *.gz > md5check.txt`
`faSize dmel-all-chromosome-r6.48.fasta.gz`
`faSize -detailed dmel-all-chromosome-r6.48.fasta.gz`
`zcat dmel-all-chromosome-r6.48.fasta.gz | grep -v '^>' | tr -d '\n' | wc -c > nucleotides.txt`
`less nucleotides.txt`
`zcat dmel-all-chromosome-r6.48.fasta.gz | grep -c '^>' > sequences.txt`
`less sequences.txt`
`zcat dmel-all-chromosome-r6.48.fasta.gz | grep -v '^>' | tr -d '\n' | grep -o N | wc -l > N.txt`
`less N.txt`
## Summary of Annotation File
`wget "http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/dmel-all-r6.48.gtf.gz" -O dmel-all-r6.48.gtf.gz`
`md5sum dmel-all-r6.48.gtf.gz`
`md5sum *gtf.gz > md5check.gtf.txt`
`less md5check.gtf.txt`
`bioawk -c gff '{print $3}' dmel-all-r6.48.gtf.gz | sort | uniq -c | sort -nr > sorted.features.txt`
`less sorted.features.txt`
`bioawk -c gff '{if ( $1 == "X" && $3 == "gene" ) print $0}' dmel-all-r6.48.gtf.gz | sort | uniq -c | wc -l > X.genes.txt`
`less X.genes.txt`
`bioawk -c gff '{if ( $1 == "2L" && $3 == "gene" ) print $0}' dmel-all-r6.48.gtf.gz | sort | uniq -c | wc -l > 2L.genes.txt`
`less 2L.genes.txt`
`bioawk -c gff '$3 == "gene" {print $1}' dmel-all-r6.48.gtf.gz | sort | uniq -c | sort -nr > Sorted.Chromosomes.txt`
`less Sorted.Chromosomes.txt`
