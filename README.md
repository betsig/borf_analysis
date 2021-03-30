# Scripts for analysis of data for testing borf

Scripts run for each species/dataset combination within their own subfolder:

human

arabidopsis

zebrafish

yeast (1 sample - main)

yeast_2sample (116M reads)

yeast_3sample (171M reads)

All Rscripts are located in the **scripts/** folder.

## Create alternative annotations for human (protein_coding only and CDS only)
Human alternative reference sets use the same Trinity assembly, but different reference .GTF/.fa files  


## All FASTQ files downloaded from EBI
```
cd fastq/
ls *_1.fastq.gz | sed 's/_1.fastq.gz//g' > filenames.txt
```

# how_are_we_stranded_here
```
cat filenames.txt | while read FILENAME; do
check_strandedness --gtf ensembl/*.gtf \
--transcripts ensembl/*.cdna.all.fa \
-r1 $FILENAME_1.fastq.gz \
-r1 $FILENAME_2.fastq.gz
done
```

# Trimgalore & Fastqc
**trimgalore-0.6.0**
**FastQC-current**
```
cat filenames.txt | while read FILENAME; do
trim_galore --fastqc --paired $FILENAME_1.fastq.gz $FILENAME_2.fastq.gz
done
```

# Rename fastq headers (Trinity likes to throw errors..)
```
cat filenames.txt | while read FILENAME; do
zcat $FILENAME_1_val_1.fq.gz | awk '{{print (NR%4 == 1) ? "@M16r1_" ++i "/1": $0}}' | gzip -c > $FILENAME_1_renamed.fq.gz
zcat $FILENAME_2_val_2.fq.gz | awk '{{print (NR%4 == 1) ? "@M16r1_" ++i "/1": $0}}' | gzip -c > $FILENAME_2_renamed.fq.gz
done
# Concatenate Fastqs
for file in *_1_renamed.fq.gz; do cat $file >> trinity_rn_1.fq.gz; done
for file in *_2_renamed.fq.gz; do cat $file >> trinity_rn_2.fq.gz; done
```

# Trinity
**Trinity-2.8.4**
**bowtie-2.3.2**
```
cd ../
Trinity --SS_lib_type RF --seqType fq --max_memory 200G --CPU 50 --output trinity_denovo --verbose --left fastq/trinity_rn_1.fq.gz --right fastq/trinity_rn_2.fq.gz
```

# kallisto
**kallisto-0.43.1**
**Trinity-2.8.4**
```
kallisto index trinity_denovo/Trinity.fasta -i trinity_kallisto_index

cat fastq/filenames.txt | while read FILENAME; do
$TRINITY_DIR/util/align_and_estimate_abundance.pl \
--seqType fq \
--left fastq/trimmed/$FILENAME_1_val_1.fq.gz \
--right fastq/trimmed/$FILENAME_2_val_2.fq.gz \
--transcripts trinity_denovo/Trinity.fasta  \
--est_method kallisto \
--trinity_mode \
--prep_reference \
--output_dir kallisto/"${FILENAME}".kallisto
done
```

## Combine kallisto counts
**kallisto_join_counts.R**

# Create Alternative Human references (CDS only/PC only)
**bedtools-2.30.0**
```
mkdir data/human_pc
mkdir data/human_cds
```

**make_alt_annotations_human.R**
Makes: pc_tx_only_ids.txt, pc_pep_only_ids.txt, human_pc.gtf, human_cds.gtf
```
gtf2bed --gtf human_cds.gtf --bed human_cds.bed
bedtools getfasta -fi Homo_sapiens.GRCh38.dna.toplevel.fa -bed human_cds.bed -name -split -s -fo human_cds.cdna.fa
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' pc_tx_only_ids.txt Homo_sapiens.GRCh38.cdna.all.fa > human_pc.cdna.fa
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' pc_pep_only_ids.txt Homo_sapiens.GRCh38.pep.all.fa > human_pc.pep.fa
cp human_pc.pep.fa human_cds.pep.fa
mv human_pc* > data/human_pc/
mv human_cds* > data/human_cds/
```

# BlastN / BlastP
**blast-2.10.1**
```
makeblastdb -in trinity_denovo/Trinity.fasta -dbtype nucl
blastn -query ../ensembl/*.cdna.all.fa -db trinity_denovo/Trinity.fasta -evalue 1e-6 -outfmt 6 > cdna2trinity_blast.txt

makeblastdb -in ensembl/*.pep1.all.fa -dbtype prot
blastx -query trinity_denovo/Trinity.fasta -db ensembl/*.pep.fa -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe qframe" > blastx_trinity2pep.txt
```

# RNAsamba
**rnasamba-0.2.4**
```
rnasamba classify rnasamba_partial.tsv trinity_denovo/Trinity.fasta ~/rnasamba/partial_length_weights.hdf5
```

# Borf of trinity
```
borf trinity_denovo/Trinity.fasta
```

# Transdecoder of trinity
**Transdecoder-5.5.0**
```
TransDecoder.LongOrfs -t trinity_denovo/Trinity.fasta
```

# Timing of Borf and Transdecoder
**make_time_tests.R**
Makes: time_orfs.sh
```
mkdir -p timing_tests
cd timing_tests
bash ./time_orfs.sh
```

**processing_timing.R**
Makes: timing_tests.txt


# Rscripts
## Functions called by other scripts:
**helper_functions.R**

## Join kallisto count files into single tsv
**kallisto_join_counts.R**
Makes: kallisto/kallisto_est_counts.tsv

## Format Ensembl/Reference annotations to have similar layouts and save as .Rdata objects
**make_ensembl_objects.R**
Makes: processed_ensembl.Rdata

## Annotate Trinity assemblies (multi-step)
**process_upto_lcs.R**
Requires: Trinity.fasta, kallisto read counts, BlastN output, BlastP output,
Reads in Trinity fasta, starts annotation & matching to best Blast hits
Makes: upto_lcs.rdata (small dataset for LCS computation), ALLDATA_upto_lcs.rdata (all objects in environment)

### Find LCS
**do_lcs_in_parts.**
Find the longest common subsequence between a translated frame and potential matched reference ORF
Run for each species dataset and frame (+1,+2,+3, -1,-2,-3) individually.
e.g.
```
cd human/
Rscript ./do_lcs_in_parts.R p1
```
This can take a really long time... not optimised code.

**recombine_lcs_data.R**
Combine the individual data from above to find LCS for each trinity/reference pair

**process_after_lcs.R**
Add the LCS data back into annotation. Finds best reference match, classifies based on ORF coverage.
Annotates ORF location (start/stop sites), UTRs, comparison to reference ORF location
Makes: trinity_translated.txt (ORF prediction on sense strand only), trinity_translated_both_strands.txt (ORF prediction on both strands -- used in manuscript)

## Annotate Trinity assemblies of alternative datasets (multi-step)
**Same process as above**, but for the human alternative reference annotations (CDS/PC only) and the different yeast assemblies (2 sample/ 3 sample)

**altdata_process_upto_lcs.R**
Reads in Trinity fasta, starts annotation & matching to best Blast hits
Makes: upto_lcs.rdata (small dataset for LCS computation), ALLDATA_upto_lcs.rdata (all objects in environment)

### Find LCS
**do_lcs_in_parts.R**
Find the longest common subsequence between a translated frame and potential matched reference ORF
Run for each species dataset and frame (+1,+2,+3, -1,-2,-3) individually.
e.g.
```
cd human_cds/
Rscript ./do_lcs_in_parts.R p1
```

This can take a really long time... not optimised code.
**altdata_recombine_lcs_data.R**
Combine the individual data from above to find LCS for each trinity/reference pair

**altdata_process_after_lcs.R**
Add the LCS data back into annotation. Finds best reference match, classifies based on ORF coverage.
Annotates ORF location (start/stop sites), UTRs, comparison to reference ORF location
Makes: trinity_translated.txt (ORF prediction on sense strand only), trinity_translated_both_strands.txt (ORF prediction on both strands -- used in manuscript)

## More data processing
**ensembl_gene_coverage.R**
Get data for coverage of the reference transcriptome by the assembled data

**ensembl_upstream_stops.R**
Makes: ensembl_upstream_stops.txt

## Make figures/tables
**assembler_citations_plot.R**
Makes: Supp_Figure_citations.pdf

**make_figures.R**
Makes: Figure1.pdf, Figure3.pdf, Figure3.pdf, Figure5.pdf,
Supp_Figure_5primes.pdf, Supp_Figure_cutoffs.pdf, Supp_Figure_read_counts.pdf

**stranded_plots.R**
Makes: Figure2.pdf

**plot_alternative_datasets.R**
Makes: Supp_Figure_altcds_altdepth.pdf, Supp_Figure_negativereplacement.pdf, Supp_Figure_genedists.pdf
