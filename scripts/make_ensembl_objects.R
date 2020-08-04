## load in GTF / GFF files

# creates: 
# ensembl_peps : peptides dataset (for lcs)
# ensembl_gff : GTF exons
# ensembl_gff_cds : gtf CDS
# ensembl_gff_gene : gtf genes (or max start/end of all gene isoforms)
# ensembl_tx_fullbody : same as above?
# utr_annot


library(data.table)
library(tidyverse)
library(stringr)
devtools::install_github("betsig/halpme")
library(halpme)
library(GenomicRanges)
library(stringdist)
library(Biostrings)
library(progress)

options(stringsAsFactors = F)
source("scripts/helper_functions.R")

for(org in c("human", "arabidopsis", "yeast","zebrafish")){
  ##### Read in 'true' annotations ######
  
  # read in ensembl annotated GFF file
  gff_file = list.files(paste0("data/reference_annotations/", org, "/"), full.names = TRUE)
  gff_file = gff_file[tolower(strv_split(gff_file, "[.]", -1)) %in% c("gff", "gff3", "gtf")][1]
  gtf_type = str_sub(toupper(strv_split(gff_file, "[.]", -1)), 1,3)
  if(org == "arabidopsis"){
    gff_file = gff_file[grep("Araport", gff_file)]
  }else{
    gff_file = gff_file[1]
  }
  
  ensembl_gff = fread(paste0("grep -v '#' ", gff_file), data.table = F, fill=T, sep='\t')
  colnames(ensembl_gff) = c('seqid','source','type','start','end','score','strand','phase','attributes')
  
  # convert roman numbers to numeric (easier for working with)
  if(org == "yeast"){
    ensembl_gff$seqid[ensembl_gff$seqid != "Mito"] = utils:::.roman2numeric(ensembl_gff$seqid[ensembl_gff$seqid != "Mito"])
  }
  
  ensembl_gff = ensembl_gff[ensembl_gff$seqid %in% c(1:99, paste0("chr", 1:99),paste0("Chr", 1:99),"ChrM", "M", "Mito", "X", "Y", "MT"),]
  
  # split into cds/gene/exon
  ensembl_gff_cds = ensembl_gff[ensembl_gff$type == "CDS",]
  if(gtf_type == "GFF" & !(org == "arabidopsis")){
    ensembl_gff_cds$transcript_id = strv_split2(ensembl_gff_cds$attributes, "transcript:", ";")
    # gene id/name later
  }else if(gtf_type == "GFF" & (org == "arabidopsis")){
    ensembl_gff_cds$transcript_id = strv_split2(ensembl_gff_cds$attributes, "Parent[=]", "[;]")
  }else if(gtf_type == "GTF"){
    ensembl_gff_cds$transcript_id = strv_split2(ensembl_gff_cds$attributes, 'transcript_id "', '";')
    ensembl_gff_cds$gene_id = strv_split2(ensembl_gff_cds$attributes, 'gene_id "', '";')
    if(all(grepl("gene_name", ensembl_gff_cds$attributes))){
      ensembl_gff_cds$gene_name = strv_split2(ensembl_gff_cds$attributes, 'gene_name "', '";')
    }else{
      ensembl_gff_cds$gene_name = ensembl_gff_cds$gene_id
    }
  }
  
  ensembl_gff_cds$width = (ensembl_gff_cds$end - ensembl_gff_cds$start) +1
  
  ensembl_gff_gene = ensembl_gff[ensembl_gff$type == "gene",]
  if(gtf_type == "GFF" & !(org == "arabidopsis")){
    ensembl_gff_gene$gene_id = strv_split2(ensembl_gff_gene$attributes , "gene[:]", ";")
  }else if(gtf_type == "GFF" & (org == "arabidopsis")){
    ensembl_gff_gene$gene_id = strv_split2(ensembl_gff_gene$attributes , "ID[=]", ";")
  }else if(gtf_type == "GTF"){
    ensembl_gff_gene$gene_id = strv_split2(ensembl_gff_gene$attributes , 'gene_id "', '";')
    if(all(grepl("gene_name", ensembl_gff_gene$attributes))){
      ensembl_gff_gene$gene_name = strv_split2(ensembl_gff_gene$attributes , 'gene_name "', '";')
    }else{
      ensembl_gff_gene$gene_name = ensembl_gff_gene$gene_id
    }
  }
  
  if(org != "yeast"){
    utr_annot = ensembl_gff[tolower(ensembl_gff$type) == "five_prime_utr" | tolower(ensembl_gff$type) == "three_prime_utr",]
  }else{
    
    ## UTR annotations (yeast only) ##
    utr_annot = fread("data/reference_annotations/yeast/GSE39128_tsedall.txt", data.table = F)
    
    # filter for utrs with >10 'counts' in at least 1 condition
    utr_annot = utr_annot[rowSums(utr_annot[,c(5:10)] > 10) > 1,]
    # make genomic ranges objects so we can filter by overlap type
    ens_granges = GRanges(seqnames = ensembl_gff_gene$seqid, strand = ensembl_gff_gene$strand, 
                          ranges = IRanges(start = ensembl_gff_gene$start, end = ensembl_gff_gene$end),
                          transcript = ensembl_gff_gene$gene_id)
    
    utr_granges = GRanges(seqnames = utr_annot$chr, strand = utr_annot$strand, 
                          ranges = IRanges(start = pmin(utr_annot$t3, utr_annot$t5), end = pmax(utr_annot$t3, utr_annot$t5)))
    
    # UTR ranges that are covering a full gene CDS region (ensembl yeast genes don't have UTRs, so we can use 'gene' annotations)
    # UTRs from this annotation cover the 5' transcription start site and the 3' transcritpion end point (i.e. over the whole gene body)
    utr_ol = findOverlaps(ens_granges, utr_granges, type = "within")
    utr_annot = get_utr_coords_from_overlaps(utr_ol, utr_granges, ens_granges)
    rm(utr_ol, utr_granges, ens_granges)
    
  }
  if(org =="arabidopsis"){
    utr_annot$transcript = strv_split2(utr_annot$attributes, "Parent[=]", ";")
    utr_annot = separate_rows(utr_annot, transcript, sep=',')
    
  }
  if(org %in% c("human", "zebrafish")){
    if(gtf_type == "GFF"){
      utr_annot$transcript = strv_split(utr_annot$attributes, "[:]",2)
    }else{
      utr_annot$transcript = remove_quotes(strv_split2(utr_annot$attributes, "transcript_id ", "[ ]"))
    }
  }
  
  utr_annot$length = abs(utr_annot$end - utr_annot$start)
  utr_annot = arrange(utr_annot, desc(length))
  utr_annot = utr_annot[!duplicated(utr_annot$transcript),]
  
  if(gtf_type == "GFF" & !(org == "arabidopsis")){
    transcripts = ensembl_gff[grepl("ID=transcript:", ensembl_gff$attributes) & 
      grepl("Parent=gene:", ensembl_gff$attributes),]
    transcripts$gene_id = strv_split2(transcripts$attributes, "gene:", ";")
    transcripts$transcript_id = strv_split2(transcripts$attributes, "transcript:", ";")
  }else if(gtf_type == "GFF" & (org == "arabidopsis")){
    transcripts = ensembl_gff[grepl("ID=", ensembl_gff$attributes) & 
                                grepl("Parent=", ensembl_gff$attributes),]
    transcripts$gene_id = strv_split2(transcripts$attributes, "Parent[=]", ";")
    transcripts$transcript_id = strv_split2(transcripts$attributes, "ID[=]", ";")
  }else if(gtf_type == "GTF"){
    transcripts = ensembl_gff[ensembl_gff$type=="transcript",]
    transcripts$gene_id = strv_split2(transcripts$attributes, 'gene_id "', '";')
    transcripts$transcript_id = strv_split2(transcripts$attributes, 'transcript_id "', '";')
  }
  
  
  ensembl_gff = ensembl_gff[ensembl_gff$type == "exon",]
  if(gtf_type == "GFF" & !(org == "arabidopsis")){
    ensembl_gff$gene_id = NA
    ensembl_gff$transcript_id = strv_split2(ensembl_gff$attributes , "transcript[:]", ";")
    # gene id/name later
  }else if(gtf_type == "GFF" & (org == "arabidopsis")){
    ensembl_gff$gene_id = NA
    ensembl_gff$transcript_id = strv_split2(ensembl_gff$attributes , "Parent[=]", ";")
    
  }else if(gtf_type == "GTF"){
    ensembl_gff$gene_id = strv_split2(ensembl_gff$attributes , 'gene_id "', '";')
    ensembl_gff$transcript_id = strv_split2(ensembl_gff$attributes , 'transcript_id "', '";')
  }
  
  if(gtf_type == "GFF"){
    ensembl_gff$gene_id = transcripts$gene_id[match(ensembl_gff$transcript_id, transcripts$transcript_id)]
    ensembl_gff_cds$gene_id = transcripts$gene_id[match(ensembl_gff_cds$transcript_id, transcripts$transcript_id)]
  }
  
  
  ensembl_tx_fullbody = aggregate(start ~ transcript_id, ensembl_gff, min)
  ensembl_tx_fullbody = cbind(ensembl_tx_fullbody, end = aggregate(start ~ transcript_id, ensembl_gff, min)[,2])
  
  ensembl_tx_fullbody = left_join(ensembl_tx_fullbody, ensembl_gff[,c(1,2,3,6,7,8,10,11)], by='transcript_id')
  ensembl_tx_fullbody = ensembl_tx_fullbody[,c(4,5,6,2,3,7,8,9,10,1)]
  ensembl_tx_fullbody = ensembl_tx_fullbody[!duplicated(ensembl_tx_fullbody$transcript_id),]
  
  # read in 'true' ensembl CDSs
  pep_file = list.files(paste0("data/reference_annotations/", org, "/"), full.names = TRUE)
  if(org == "arabidopsis"){
    pep_file = pep_file[grepl("pep1", pep_file) & grepl("Araport", pep_file)]
    pep_file = pep_file[str_sub(pep_file, -10,-1) == "pep1.fasta"]
  }else{
    pep_file = pep_file[str_sub(pep_file, -11,-1) == "pep1.all.fa"]
  }
  pep_file = pep_file[!grepl("blast", pep_file)]
  
  ensembl_peps = halpme::read_fasta2df(pep_file)
  if(org == "arabidopsis"){
    ensembl_peps$transcript = strv_split(ensembl_peps$seq_id, "[ ][|]", 1)
  }else{
    ensembl_peps$transcript = strv_split2(ensembl_peps$seq_id, "transcript:", "[ ]")
  }
  
  ensembl_peps$length = nchar(ensembl_peps$seq)
  
  save(ensembl_peps, ensembl_tx_fullbody, ensembl_gff, ensembl_gff_cds,ensembl_gff_gene, transcripts, utr_annot, file=paste0("data/", org, "/processed_ensembl.Rdata"))
}
