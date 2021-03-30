## Clean and ANNOTATED/FUNCTIONS for yeast borf comps

library(data.table)
library(tidyverse)
library(stringr)
library(halpme)
library(GenomicRanges)
library(stringdist)
library(Biostrings)

options(stringsAsFactors = F)
source("scripts/helper_functions.R")


alt_orgs = c('human', 'human', 'yeast', 'yeast')
alt_versions = c("cds", "pc", "2sample", "3sample")

for(i in 1:length(alt_orgs)){
  org = alt_orgs[i]
  version = alt_versions[i]
  org_version = ifelse(version=="", org, paste(org, version, sep="_"))
  org_version_assembly_data = ifelse(org == "human", org, org_version)
  
  if(!file.exists(paste0("data/",org_version, "/ALLDATA_upto_lcs.rdata"))){
  
    load(paste0("data/", org, "/processed_ensembl.Rdata"))
      
    #### read in Trinity assembly fasta ####
    trinity = halpme::read_fasta2df(paste0("data/", org_version_assembly_data, "/Trinity.fasta"))
    trinity$transcript_id = strv_split(trinity$seq_id, "[ ]", 1)
    trinity$gene_id = strv_split(trinity$transcript_id, "_i", 1)
    trinity$seq_len = nchar(trinity$seq)
    #trinity$seq = NULL # keep seq for LCS stuff later

    # kallisto read counts
    make_kallisto_summary(org_version = org_version)
    kallisto = fread(paste0("data/", org_version_assembly_data, "/kallisto/kallisto_est_counts.tsv"),data.table = F, fill=T, sep='\t')
    kallisto$total_counts = rowSums(kallisto[,-1])
    
    combined_data = trinity[,-which(colnames(trinity) =='seq')]
    combined_data$counts= kallisto$total_counts[match(combined_data$transcript_id, kallisto$target_id)]
    combined_data$reads_per_kb = (combined_data$counts / combined_data$seq_len) * 1000
    
    write.table(combined_data, file=paste0("data/", org_version, "/trinity_with_kallisto.txt"), sep='\t', quote=F, row.names = F)
  
    # blast cdna --> trinity
    e2t = fread(paste0("data/",org_version, "/blast_results/cdna2trinity_blast.txt"), data.table=F)
    colnames(e2t) = c("query", "subject", "p_ident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    
    # check blastx run on protein ids???
    ensembl_peps$protein_id = strv_split(ensembl_peps$seq_id, "[ ]", 1)
    ensembl_peps$first_aa = ifelse(str_sub(ensembl_peps$seq, 1,1) == "M", "M", "non-M")
    
    ###########################
    
    e2t = arrange(e2t, evalue)
    # filter for 1e-30
    e2t = e2t[e2t$evalue < 1e-30,]
    e2t$trinity_gene = strv_split(e2t$subject, "_i",1)
    e2t$read_counts = combined_data$counts[match(e2t$subject, combined_data$transcript_id)]
    e2t$reads_per_kb = combined_data$reads_per_kb[match(e2t$subject, combined_data$transcript_id)]
    max_rpk = aggregate(reads_per_kb ~ gene_id, combined_data, max)
    e2t$gene_rpk = max_rpk$reads_per_kb[match(e2t$trinity_gene, max_rpk$gene_id)]
    e2t$match_strand = ifelse(((e2t$qend - e2t$qstart > 0) & (e2t$send - e2t$sstart > 0)) | ((e2t$qend - e2t$qstart < 0) & (e2t$send - e2t$sstart < 0)), "+", "-")
    e2t$query_gene = gsub("_mRNA", "", e2t$query)
    rm(max_rpk)
    
    # filter by read counts (min 100 reads per 1000 nt of assembled 'transcript' or an average of 10x coverage for 100bp reads)
    e2t.filtered = arrange(e2t, evalue)
    e2t.filtered = e2t.filtered[e2t.filtered$reads_per_kb >= 100 | e2t.filtered$read_counts >=10,]
    e2t.filtered = e2t.filtered[!duplicated(paste0(e2t.filtered$query, e2t.filtered$subject)),]
    # split into trinity genes (not isoforms) with a single blast hit and those with multiple hits
    no_dups = which(!duplicated(paste0(e2t.filtered$query, e2t.filtered$trinity_gene)))
    one_match_trinity_genes = table(e2t.filtered$trinity_gene[no_dups]) %>% as.data.frame() %>% filter(Freq==1)
    e2t_single_gene_match = e2t.filtered[e2t.filtered$trinity_gene %in% one_match_trinity_genes$Var1,]
    e2t_multi_gene_match = e2t.filtered[!(e2t.filtered$trinity_gene %in% one_match_trinity_genes$Var1),]
    rm(e2t.filtered, no_dups, one_match_trinity_genes)
    
    # POSITIVE STRAND matches ONLY
    # filter by read counts (min 100 reads per 1000 nt of assembled 'transcript' or an average of 10x coverage for 100bp reads)
    e2t.filtered = arrange(e2t, evalue)
    e2t.filtered = e2t.filtered[e2t.filtered$reads_per_kb >= 100 | e2t.filtered$read_counts >=10,]
    e2t.filtered = e2t.filtered[e2t.filtered$match_strand == "+",]
    e2t.filtered = e2t.filtered[!duplicated(paste0(e2t.filtered$query, e2t.filtered$subject)),]
    # split into trinity genes (not isoforms) with a single blast hit and those with multiple hits
    no_dups = which(!duplicated(paste0(e2t.filtered$query, e2t.filtered$trinity_gene)))
    one_match_trinity_genes = table(e2t.filtered$trinity_gene[no_dups]) %>% as.data.frame() %>% filter(Freq==1)
    e2t_single_gene_match_pos = e2t.filtered[e2t.filtered$trinity_gene %in% one_match_trinity_genes$Var1,]
    e2t_multi_gene_match_pos = e2t.filtered[!(e2t.filtered$trinity_gene %in% one_match_trinity_genes$Var1),]
    
    # write trinity ids to file so we don't blast transcripts that aren't expressed
    write.table(unique(c(e2t_multi_gene_match$subject, e2t_multi_gene_match_pos$subject)), 
                file=paste0("data/", org_version, "/blastx_filter_trinity_ids.txt"), 
                quote=F, row.names = F, col.names = F, sep='\n')
    
    # blastx to match to best
    blastx = fread(paste0("data/",org_version, "/blast_results/blastx_trinity2pep.txt"), data.table = F, fill=T)
    colnames(blastx) = c("query", "subject", "p_ident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","sframe", "qframe")
    
    ## combine multi blastx hits to single line
    blastx = combine_blastx_lines(blastx)
    
    blastx$transcript_id = ensembl_peps$transcript[match(blastx$subject, ensembl_peps$protein_id)]
    blastx$ensembl_pep_len = ensembl_peps$length[match(blastx$transcript_id, ensembl_peps$transcript)]
    blastx$subject_match_length = abs(blastx$send - blastx$sstart) +1
    blastx$subject_first_aa = ensembl_peps$first_aa[match(blastx$transcript_id, ensembl_peps$transcript)]
    blastx$subject_coverage = blastx$subject_match_length / blastx$ensembl_pep_len
    
    blastx = arrange(blastx, query, evalue, desc(subject_coverage), desc(p_ident*length))
    blastx$match_strand = ifelse(blastx$qframe < 0, "-","+")
    
    
    e2t_multi_gene_match = e2t_multi_gene_match[,c(1:18)]
    e2t_multi_gene_match = left_join(e2t_multi_gene_match, blastx, suffix = c('', '.blastx'), by=c('subject'='query', 'query'='transcript_id'))
    
    # filter by subject coverage FIRST, then evalue etc...
    e2t_multi_gene_match = arrange(e2t_multi_gene_match, subject, desc(subject_coverage), evalue.blastx, 
                                   desc(p_ident.blastx*length.blastx), evalue, desc(bitscore))
    multi_match_best_hit = e2t_multi_gene_match[!duplicated(e2t_multi_gene_match$subject),]
    multi_match_best_hit$match_type = "multi"
    e2t_single_gene_match$match_type = "single"
    
    e2t_best_matches = rbind(e2t_single_gene_match[,c(1,2)], multi_match_best_hit[,c(1,2)])
    
    e2t_multi_gene_match_pos = e2t_multi_gene_match_pos[,c(1:18)]
    e2t_multi_gene_match_pos = left_join(e2t_multi_gene_match_pos, blastx, suffix = c('', '.blastx'), by=c('subject'='query', 'query'='transcript_id'))
    
    e2t_multi_gene_match_pos = arrange(e2t_multi_gene_match_pos, subject, desc(subject_coverage),evalue.blastx, 
                                       desc(p_ident.blastx*length.blastx), evalue, desc(bitscore))
    multi_match_best_hit_pos = e2t_multi_gene_match_pos[!duplicated(e2t_multi_gene_match_pos$subject),]
    multi_match_best_hit_pos$match_type = "multi"
    e2t_single_gene_match_pos$match_type = "single"
    e2t_best_matches_pos = rbind(e2t_single_gene_match_pos[,c(1,2)], multi_match_best_hit_pos[,c(1,2)])
    
    e2t_best_matches_do_pos = rbind(e2t_best_matches, e2t_best_matches_pos) %>% distinct()
    e2t_best_matches_do_neg = e2t_best_matches
    
    save(e2t_best_matches_do_pos, e2t_best_matches_do_neg, trinity, ensembl_peps,longest_common_substring, file=paste0("data/", org_version, "/upto_lcs.rdata"))
    save.image(paste0("data/",org_version, "/ALLDATA_upto_lcs.rdata"))
  }
}

#####################################################################
# then run lcs in parts