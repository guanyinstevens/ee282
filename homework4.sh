## Commands for Homework 4
# Calculate the following for all sequences ≤ 100kb and all sequences > 100kb:
`wget "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/fasta/dmel-all-chromosome-r6.60.fasta.gz" -O dmel-all-chromosome-r6.60.fasta.gz`
`wget "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/fasta/md5sum.txt" -O md5sum.txt`
`md5sum -c md5sum.txt`
`faFilter minSize=100000 dmel-all-chromosome-r6.60.fasta.gz dmel-all-chromosome-r6.60_more_100kb.fasta`
`faSize dmel-all-chromosome-r6.60_more_100kb.fasta`
`faFilter maxSize=100000 dmel-all-chromosome-r6.60.fasta.gz dmel-all-chromosome-r6.60_less_100kb.fasta`
`faSize dmel-all-chromosome-r6.60_less_100kb.fasta`
`bioawk -c fastx '{ print $name, length($seq) }' dmel-all-chromosome-r6.60_less_100kb.fasta > length-lte-100kb-dmel-all-chromosome-r6.60.txt`
`bioawk -c fastx '{ print $name, length($seq) }' dmel-all-chromosome-r6.60_more_100kb.fasta > length-gt-100kb-dmel-all-chromosome-r6.60.txt`
`bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.60_less_100kb.fasta > gc_lte_100kb_dmel-all-chromosome-r6.60.txt`
`bioawk -c fastx '{ print $name,gc($seq) }' dmel-all-chromosome-r6.60_more_100kb.fasta > gc_gt_100kb_dmel-all-chromosome-r6.60.txt`
# Plots of the following for for all sequences ≤ 100kb and all sequences > 100kb:
In R, unless otherwise specified. 
`gc_content_lte_100kb.df <- read.delim('gc_lte_100kb_dmel-all-chromosome-r6.60.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=NULL)`
`gc_content_lte_100kb.df <- gc_content_lte_100kb.df[, c("Scaffold", "X")]`
`"ggplot(gc_content_lte_100kb.df, aes(x = Scaffold, y = X)) +
geom_bar(stat = ""identity"", color = ""black"") + labs(title = ""GC Content by Scaffold"", x = ""Scaffold"", y = ""GC Content (%)"")"`
`gc_content_gte_100kb.df <- read.delim('gc_gt_100kb_dmel-all-chromosome-r6.60.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names=NULL)`
`ggplot(gc_content_gte_100kb.df, aes(x = V1, y = V2)) + geom_bar(stat = ""identity"")`
`ggplot(gc_content_gte_100kb.df, aes(x = V1, y = V2)) + geom_bar(stat = ""identity"", fill = ""pink"", color = ""black"") + labs(title = ""GC Content for Sequences Larger than 100 kb"", x = ""Scaffold"", y = ""GC Content (%)"")`
`length_gte_100kb.df <- read.delim('length-gt-100kb-dmel-all-chromosome-r6.60.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)`
`length_gte_100kb.df$log_V2 <- log10(length_gte_100kb.df$V2)
ggplot(length_gte_100kb.df, aes(x = V1, y = log_V2, fill = V1)) + geom_col() + labs(title = ""Log-Transformed Sequence Lengths by Scaffold"", x = ""Scaffold"", y = ""Log10 Sequence Length"")`
`length_lte_100kb.df <- read.delim('length-lte-100kb-dmel-all-chromosome-r6.60.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)`
`length_lte_100kb.df$log_V2 <- log10(length_lte_100kb.df$V2)`
`ggplot(length_lte_100kb.df, aes(x = log_V2)) + geom_histogram(bins = 50, fill = "lavender", color = "black") +  labs(title = "Histogram of Sequence Lengths Less than 100kb", x = "Log10 Sequence Length", y = "Frequency of Sequence Length")`
Using terminal
`#!/usr/bin/env bash`
`sort -rnk2 length-lte-100kb-dmel-all-chromosome-r6.60.txt > sorted_length-lte-100kb-dmel-all-chromosome-r6.60.txt`
`sort -rnk2 length-gt-100kb-dmel-all-chromosome-r6.60.txt > sorted_length-gt-100kb-dmel-all-chromosome-r6.60.txt`
`plotCDF sorted_length-gt-100kb-dmel-all-chromosome-r6.60.txt sorted_length_gt_100kb_CDF.png`
`plotCDF sorted_length-lte-100kb-dmel-all-chromosome-r6.60.txt sorted_length_lte_100kb_CDF.png`

# Genome Assembly
`hifiasm -o ISO1.asm -t16 -l0 ISO1_Hifi_AdaptorRem.40X.fasta.gz 2> ISO1.asm.log`
`awk '$1 == "S" {print ">" $2 "\n" $3}' ISO1.asm.bp.p_ctg.gfa > ISO1.asm.bp.p_ctg.fasta`
# N50 Calculations
 `bioawk -c fastx '{print $name, length($seq)}' ISO1.asm.bp.p_ctg.fasta | sort -n -k2,2 -r > sorted_sequences_by_length.txt`
`total_length=$(awk '{sum += $2} END {print sum}' sorted_sequences_by_length.txt)`
`half_total_length=$(echo "$total_length / 2" | bc)`
`awk -v half_total_length="$half_total_length" ' BEGIN {cumulative_length = 0}{ cumulative_length += $2; if (cumulative_length >= half_total_length && n50 == 0) {n50 = $2; print "N50: " n50; exit} }' sorted_sequences_by_length.txt`
`wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz" -O GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz`
`wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/md5checksums.txt" -O md5checksums.txt`
`md5sum -c md5checksums.txt`
`plotCDF dmel_contig_sorted_sequences_by_length.txt dmel_sorted_sequences_by_length.txt sorted_sequences_by_length.txt assembly_comparison.png`
`zcat GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | grep -v ">" | tr -d "\n" | grep -o "N*" | awk '{ print length }' > n_lengths.txt`
`wc -l n_lengths.txt`
`aSplitByN GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz contig_GCF_Release_6_ISO1_genomic.fna.gz 572`
`bioawk -c fastx '{print $name, length($seq)}' contig_GCF_Release_6_ISO1_genomic.fna.gz | sort -n -k2,2 -r > contig_sorted_sequences_by_length.txt `
`contig_total_length=$(awk '{sum += $2} END {print sum}' contig_sorted_sequences_by_length.txt)`
`contig_half_total_length=$(echo "$contig_total_length / 2" | bc)`
`awk -v contig_half_total_length=$contig_half_total_length"" ' BEGIN {cumulative_length = 0}{ cumulative_length += $2; if (cumulative_length >= contig_half_total_length && n50 == 0) {n50 = $2; print "N50: " n50; exit} }' contig_sorted_sequences_by_length.txt`
`bioawk -c fastx '{print $name, length($seq)}' GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | sort -n -k2,2 -r > scaffold_sorted_sequences_by_length.txt`
`scaffold_total_length=$(awk '{sum += $2} END {print sum}' scaffold_sorted_sequences_by_length.txt)`
`scaffold_half_total_length=$(echo "$scaffold_total_length / 2" | bc)`
`awk -v scaffold_half_total_length="$scaffold_half_total_length" ' BEGIN {cumulative_length = 0}{ cumulative_length += $2; if (cumulative_length >= scaffold_half_total_length && n50 == 0) {n50 = $2; print "N50: " n50; exit} }' scaffold_sorted_sequences_by_length.txt`
`plotCDF contig_sorted_sequences_by_length.txt scaffold_sorted_sequences_by_length.txt sorted_sequences_by_length.txt assembly.png`
# BUSCO Scores
`mamba install compleasm`
`mamba update compleasm`
`mamba create --name compleasm_env compleasm`
`compleasm -download eukaryota_odb10`
`compleasm run -a ISO1.asm.bp.p_ctg.fasta -o /pub/glsteven -t '16' -l 'eukaryota_odb10'`
`compleasm run -a contig_GCF_Release_6_ISO1_genomic.fna.gz -o /pub/glsteven -t '16' -l 'eukaryota_odb10'`
`compleasm run -a GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz -o /pub/glsteven -t '16' -l 'eukaryota_odb10'`
`compleasm download diptera_odb10`
`compleasm run -a ISO1.asm.bp.p_ctg.fasta -o /pub/glsteven -t '16' -l 'diptera_odb10'`
`compleasm run -a contig_GCF_Release_6_ISO1_genomic.fna.gz -o /pub/glsteven -t '16' -l 'diptera_odb10'`
`compleasm run -a GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz -o /pub/glsteven -t '16' -l 'diptera_odb10'`

