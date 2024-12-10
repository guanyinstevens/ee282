# Homework 4

## Downloading data and performing checksums
To download the data, the following commands and checks were performed. 
`wget "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/fasta/dmel-all-chromosome-r6.60.fasta.gz" -O dmel-all-chromosome-r6.60.fasta.gz`
`wget "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/fasta/md5sum.txt" -O md5sum.txt`
`md5sum -c md5sum.txt`
The checksum should fail for all other fiels, but say okay for the dmel-all-chromosome-r6.60.fasta.gz
## Calculating the number of nucleotides, sequences, and N's that are ≤ 100kb and > 100kb (less than or equal to & greater than 100kb)
`faFilter minSize=100000 dmel-all-chromosome-r6.60.fasta.gz dmel-all-chromosome-r6.60_more_100kb.fasta`
This command will sort the input file (our fasta file that contains sequences for all chromosomes) by files that are greater than 100kb. 
`faSize dmel-all-chromosome-r6.60_more_100kb.fasta`
Performing faSize on the sorted file will tell us the number of bases/nucleotides, Ns, and sequences. 
The number of bases that are above 100kb is 137547960. 
The number of Ns that are above 100kb is 490385. 
The total number of sequences above 100kb is 7.  
`faFilter maxSize=100000 dmel-all-chromosome-r6.60.fasta.gz dmel-all-chromosome-r6.60_less_100kb.fasta`
This will sort the input file of all the chromosomes in the genome that are less than or equal to 100kb. 
`faSize dmel-all-chromosome-r6.60_less_100kb.fasta`
This will give an output with details of the sequences that are less than or equal to 100kb. From this we can look at the number of bases/nucleotides, Ns, and sequences. 
The number of bases that are less than or equal to 100kb: 6178042
The number of Ns that are less than or equal to 100kb: 662593
The total number of sequences that are less than or equal to 100kb: 1863. 
## Plots of the following sequences <= 100kb and all sequences > 100kb:
First, to be able to plot everything, I created txt files from the original fasta files using the following commands:
`bioawk -c fastx '{ print $name, length($seq) }' dmel-all-chromosome-r6.60_less_100kb.fasta > length-lte-100kb-dmel-all-chromosome-r6.60.txt`
`bioawk -c fastx '{ print $name, length($seq) }' dmel-all-chromosome-r6.60_more_100kb.fasta > length-gt-100kb-dmel-all-chromosome-r6.60.txt`
`bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.60_less_100kb.fasta > gc_lte_100kb_dmel-all-chromosome-r6.60.txt`
`bioawk -c fastx '{ print $name,gc($seq) }' dmel-all-chromosome-r6.60_more_100kb.fasta > gc_gt_100kb_dmel-all-chromosome-r6.60.txt`
I then moved those files to my home computer to be able to make the plots using R studio. 
To download the data into R, I used the command: `gc_content_lte_100kb.df <- read.delim('gc_lte_100kb_dmel-all-chromosome-r6.60.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=NULL)`
I then sorted the columns for the scaffolds and GC content (which were labeled as X) by using the command: `> gc_content_lte_100kb.df <- gc_content_lte_100kb.df[, c("Scaffold", "X")]`. 
To then plot the data I used the command: `"ggplot(gc_content_lte_100kb.df, aes(x = Scaffold, y = X)) +
geom_bar(stat = ""identity"", color = ""black"") + labs(title = ""GC Content by Scaffold"", x = ""Scaffold"", y = ""GC Content (%)"")"`.
Plot: 
![GC Content for Less than 100kb](https://github.com/guanyinstevens/ee282/blob/homework4/GC_content_lte_100kb.png?raw=true)
[here](https://github.com/guanyinstevens/ee282/blob/homework4/GC_content_lte_100kb.png)
This plot would be the plot of the GC content for sequences that were less than or equal to 100kb. 
To download the file that contains the sequences that are greater than 100kb, I used the command: `gc_content_gte_100kb.df <- read.delim('gc_gt_100kb_dmel-all-chromosome-r6.60.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=NULL)`
To plot the GC content for the sequences that were greater than 100kb, I used the following commands:
`ggplot(gc_content_gte_100kb.df, aes(x = V1, y = V2)) + geom_bar(stat = ""identity"")`
`ggplot(gc_content_gte_100kb.df, aes(x = V1, y = V2)) + geom_bar(stat = ""identity"", fill = ""pink"", color = ""black"") + labs(title = ""GC Content for Sequences Larger than 100 kb"", x = ""Scaffold"", y = ""GC Content (%)"")`
Plot: 
![Plot of GC Content for Greater than 100kb](https://github.com/guanyinstevens/ee282/blob/homework4/Plot_of_GC_Content_for_Greater_than_100kb.png?raw=true)
[here](https://github.com/guanyinstevens/ee282/blob/homework4/Plot_of_GC_Content_for_Greater_than_100kb.png)
For sequence length, I plotted log10 values of the sequences to better visualize the frequency of the lengths. 
For sequences that were greater than 100kb, I used the following commands to plot the sequence lengths:
`length_gte_100kb.df <- read.delim('length-gt-100kb-dmel-all-chromosome-r6.60.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)`
`length_gte_100kb.df$log_V2 <- log10(length_gte_100kb.df$V2)
ggplot(length_gte_100kb.df, aes(x = V1, y = log_V2, fill = V1)) + geom_col() + labs(title = ""Log-Transformed Sequence Lengths by Scaffold"", x = ""Scaffold"", y = ""Log10 Sequence Length"")`
Plot: 
![Log Transformed Sequence Length (>=100kb)](https://github.com/guanyinstevens/ee282/blob/homework4/log_transformed_sequence_length_gte_100kb.png?raw=true)
[here](https://github.com/guanyinstevens/ee282/blob/homework4/log_transformed_sequence_length_gte_100kb.png)
To then plot sequences that were less than or equal to 100kb, I used the following commands:
`length_lte_100kb.df <- read.delim('length-lte-100kb-dmel-all-chromosome-r6.60.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)`
`length_lte_100kb.df$log_V2 <- log10(length_lte_100kb.df$V2)`
`ggplot(length_lte_100kb.df, aes(x = log_V2)) + geom_histogram(bins = 50, fill = "lavender", color = "black") +  labs(title = "Histogram of Sequence Lengths Less than 100kb", x = "Log10 Sequence Length", y = "Frequency of Sequence Length")`
Plot: 
![Log Transformed Sequence Length (<=100kb)](https://github.com/guanyinstevens/ee282/blob/homework4/log_transformed_sequence_length_lte_100kb.png?raw=true)
[here](https://github.com/guanyinstevens/ee282/blob/homework4/log_transformed_sequence_length_lte_100kb.png)
The commands then produce a histogram showing the log10 values of the different sequences that were greater than or less than or equal to 100kb. 
To plot the cumulative sequence size from largest to smallest sequences, I used the command line in Terminal, using the plotCDF utility. I first used the sort command to sort the sequence length from largest to smallest, using the following commands:
`sort -rnk2 length-lte-100kb-dmel-all-chromosome-r6.60.txt > sorted_length-lte-100kb-dmel-all-chromosome-r6.60.txt`
`sort -rnk2 length-gt-100kb-dmel-all-chromosome-r6.60.txt > sorted_length-gt-100kb-dmel-all-chromosome-r6.60.txt`
To the plot the cumulative sequence sizes, I used the commands:
`plotCDF sorted_length-gt-100kb-dmel-all-chromosome-r6.60.txt sorted_length_gt_100kb_CDF.png`
`plotCDF sorted_length-lte-100kb-dmel-all-chromosome-r6.60.txt sorted_length_lte_100kb_CDF.png`
Plot: 
![Sorted Length > 100kb CDF](https://github.com/guanyinstevens/ee282/blob/homework4/sorted_length_gt_100kb_CDF.png?raw=true)
[here](https://github.com/guanyinstevens/ee282/blob/homework4/sorted_length_gt_100kb_CDF.png)
![Sorted Length <= 100kb CDF](https://github.com/guanyinstevens/ee282/blob/homework4/sorted_length_lte_lte_100kb_CDF.png?raw=true)
[here](https://github.com/guanyinstevens/ee282/blob/homework4/sorted_length_lte_lte_100kb_CDF.png)
## Genome Assembly
To first download and perform the hifiasm assembly, I did the following commands:
`hifiasm -o ISO1.asm -t16 -l0 ISO1_Hifi_AdaptorRem.40X.fasta.gz 2> ISO1.asm.log`
To then calculate the N50 value for the assembly, I needed to sort the scaffold by length from largest to smallest. I did this by first creating a fasta file from the contig gfa file, using this command: `awk '$1 == "S" {print ">" $2 "\n" $3}' ISO1.asm.bp.p_ctg.gfa > ISO1.asm.bp.p_ctg.fasta`
To then perform the sorting, I used bioawk and sort: `bioawk -c fastx '{print $name, length($seq)}' ISO1.asm.bp.p_ctg.fasta | sort -n -k2,2 -r > sorted_sequences_by_length.txt`
For N50 calculations, I needed to calculate the total length, the median point of the total length, and then calculate the cumulative length. This was done using the following commands:
`total_length=$(awk '{sum += $2} END {print sum}' sorted_sequences_by_length.txt)`
`half_total_length=$(echo "$total_length / 2" | bc)`
`awk -v half_total_length="$half_total_length" ' BEGIN {cumulative_length = 0}{ cumulative_length += $2; if (cumulative_length >= half_total_length && n50 == 0) {n50 = $2; print "N50: " n50; exit} }' sorted_sequences_by_length.txt`
This should have an output of N50: 21715751
On the website for the Genome assembly Release 6 plus ISO1 MT, their contig N50 is 21.5 Mb. My assembly was 21.7 Mb, which would be a little larger and more contiguous. 
To the compare the N50's for the scaffold and contig sequences, I first downloaded the relevant data and performed a md5 check. 
`wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz" -O GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz`
`wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/md5checksums.txt" -O md5checksums.txt`
`md5sum -c md5checksums.txt`
To assemble just the contigs, I used faSplitByN. 
I split by 10 N's using this command:
`FaSplitByN GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz contig_GCF_Release_6_ISO1_genomic.fna.gz 20` gave me my fasta file for just the contigs. 
I then used the following commands to calculate the N50 for the scaffold assembly:
`bioawk -c fastx '{print $name, length($seq)}' contig_GCF_Release_6_ISO1_genomic.fna.gz | sort -n -k2,2 -r > contig_sorted_sequences_by_length.txt ` which allowed me to sort my sequences from largest to smallest. 
To calculate the full length and half length, I used the commands:
`contig_total_length=$(awk '{sum += $2} END {print sum}' contig_sorted_sequences_by_length.txt)`
`contig_half_total_length=$(echo "$contig_total_length / 2" | bc)`
To then calculate the N50, I ran this command: `awk -v contig_half_total_length=$contig_half_total_length"" ' BEGIN {cumulative_length = 0}{ cumulative_length += $2; if (cumulative_length >= contig_half_total_length && n50 == 0) {n50 = $2; print "N50: " n50; exit} }' contig_sorted_sequences_by_length.txt`
The N50 for the contig was: 23513712 (23.5 Mb)
For the scaffold, I ran the same commands:
`bioawk -c fastx '{print $name, length($seq)}' GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | sort -n -k2,2 -r > scaffold_sorted_sequences_by_length.txt`
`scaffold_total_length=$(awk '{sum += $2} END {print sum}' scaffold_sorted_sequences_by_length.txt)`
`scaffold_half_total_length=$(echo "$scaffold_total_length / 2" | bc)`
`awk -v scaffold_half_total_length="$scaffold_half_total_length" ' BEGIN {cumulative_length = 0}{ cumulative_length += $2; if (cumulative_length >= scaffold_half_total_length && n50 == 0) {n50 = $2; print "N50: " n50; exit} }' scaffold_sorted_sequences_by_length.txt`
The N50 score was: 25286936 (25.2 Mb)


To plot the cumulative sequences sizes, I used plotCDF and used the following command: `plotCDF contig_sorted_sequences_by_length.txt scaffold_sorted_sequences_by_length.txt sorted_sequences_by_length.txt assembly.png`
Plot: 
![Assembly Plot](https://github.com/guanyinstevens/ee282/blob/homework4/assembly.png?raw=true)
[here](https://github.com/guanyinstevens/ee282/blob/homework4/assembly.png)

## Busco Analysis
To perform the BUSCO analysis, I tried to use the BUSCO package and pipline but was unable to. I instead used the compleasm package as an alternative. 
The commands I used to install the packages were:
`mamba install compleasm`
`mamba update compleasm`
`mamba create --name compleasm_env compleasm`
I then downloaded the eukaryotic lineage to perform my analysis, using the following command: `compleasm -download eukaryota_odb10`
I first analyzed the ISO1 hifiasm assembly using the following command: `compleasm run -a ISO1.asm.bp.p_ctg.fasta -o /pub/glsteven -t '16' -l 'eukaryota_odb10'`
The results I got were as follows:
S:100.00%, 255
D:0.00%, 0
F:0.00%, 0
I:0.00%, 0
M:0.00%, 0
N:255
 
To then analyze the contig assembly, I used a similar command: `compleasm run -a contig_GCF_Release_6_ISO1_genomic.fna.gz -o /pub/glsteven -t '16' -l 'eukaryota_odb10'`
The results from the command were: 
S:100.00%, 255
D:0.00%, 0
F:0.00%, 0
I:0.00%, 0
M:0.00%, 0
N:255

To then analyze the scaffold assembly, I used the following command: `compleasm run -a GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz -o /pub/glsteven -t '16' -l 'eukaryota_odb10'`
The results were:
S:100.00%, 255
D:0.00%, 0
F:0.00%, 0
I:0.00%, 0
M:0.00%, 0
N:255

This indicates for all three assemblies, that the number of BUSCO genes that aligned once to the genome was 255, and the number of unmapped genes was 255. There were no duplciated genes or fragments genes. However, the percentage for the matched genes was 100%, so I am unsure why there was unmapped genes. Maybe it is because it is an inbred strain of *Drosophila* and not a true wild-type genome, so the core markers are present, while other genes may not be found. 

In addition, I also performed the same analysis, but instead of using the eukaryote lineage, using the diptera lineage. 
I first downloaded the diptera lineage from the remote repository using the following command: `compleasm download diptera_odb10`
I then analyzed the hifiasm assembly, contig assembly, and scaffold assembly using the respective commands: 
hifiasm - `compleasm run -a ISO1.asm.bp.p_ctg.fasta -o /pub/glsteven -t '16' -l 'diptera_odb10'`
contig - `compleasm run -a contig_GCF_Release_6_ISO1_genomic.fna.gz -o /pub/glsteven -t '16' -l 'diptera_odb10'`
scaffold - `compleasm run -a GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz -o /pub/glsteven -t '16' -l 'diptera_odb10'`
The results from each command were as follows:
Hifiasm: 
S:99.73%, 3276
D:0.27%, 9
F:0.00%, 0
I:0.00%, 0
M:0.00%, 0
N:3285

Contig:
S:99.73%, 3276
D:0.27%, 9
F:0.00%, 0
I:0.00%, 0
M:0.00%, 0
N:3285

Scaffold: 
S:99.73%, 3276
D:0.27%, 9
F:0.00%, 0
I:0.00%, 0
M:0.00%, 0
N:3285

This tells us that the number of BUSCO genes that were successfully aligned to the lineage genome was 99.73%, with 3276 genes entirely aligned once. There were 9 duplicated genes to the diptera lineage, but no fragmented genes were present. There were 3285 missing genes, which again is interesting because of the high single copy percentage. 
