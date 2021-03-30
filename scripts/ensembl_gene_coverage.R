
library(data.table)
library(tidyverse)
library(stringr)
library(halpme)
library(GenomicRanges)
library(stringdist)
library(Biostrings)

options(stringsAsFactors = F)

make_ens_coverage = function(org, ensembl_gff_cds, tt){
    ens_gene_cov = data.frame(gene_id = unique(ensembl_gff_cds$gene_id))
    ens_gene_cov$class = "not_assembled"
    ens_gene_cov$class[ens_gene_cov$gene_id %in% tt$gene_id] = "assembled"
    ens_gene_cov$class[ens_gene_cov$gene_id %in% tt$gene_id[tt$read_counts >= 10]] = "assembled_count_10"
    ens_gene_cov$class[ens_gene_cov$gene_id %in% tt$gene_id[tt$read_counts >= 10 & tt$base_class %in% c("complete", "incomplete_3prime", "incomplete", "incomplete_5prime")]] = "assembled_count_10_ORF_match"
    ens_gene_cov$class[ens_gene_cov$gene_id %in% tt$gene_id[tt$read_counts >= 100]] = "assembled_count_100"
    ens_gene_cov$class[ens_gene_cov$gene_id %in% tt$gene_id[tt$read_counts >= 100 & tt$base_class %in% c("complete", "incomplete_3prime", "incomplete", "incomplete_5prime")]] = "assembled_count_100_ORF_match"
    ens_gene_cov$org = org
    return(ens_gene_cov)
}



orgs = c('arabidopsis','yeast', 'zebrafish','human')

ens_gene_cov = NULL
for(i in seq_along(orgs)){
    org = orgs[i]
    
    tt_pos = fread(paste0("data/", org, "/trinity_translated.txt"), data.table = F)
    tt_pos$org = org
    
    load(paste0("data/", org, "/processed_ensembl.Rdata"))
    if(org == "yeast"){
        tt_pos$gene_id = ensembl_gff_cds$gene_id[match(strv_split(tt_pos$query, "[.]",1), ensembl_gff_cds$transcript_id)]
    }else if(org=="human"){
        tt_pos$gene_id = ensembl_gff_cds$gene_id[match(strv_split(tt_pos$query_gene, "[.]",1), ensembl_gff_cds$transcript_id)]
    }else if(org=="arabidopsis"){
        tt_pos$gene_id = ensembl_gff_cds$gene_id[match((tt_pos$query_gene), ensembl_gff_cds$transcript_id)]
    }else if(org=="zebrafish"){
        tt_pos$gene_id = ensembl_gff_cds$gene_id[match(strv_split(tt_pos$query_gene, "[.]",1), ensembl_gff_cds$transcript_id)]
    }
    
    ens_gene_cov_org = make_ens_coverage(org, ensembl_gff_cds, tt_pos)
    ens_gene_cov = rbind(ens_gene_cov,ens_gene_cov_org)
    
}

ens_gene_cov$class = factor(ens_gene_cov$class, levels = rev(c("assembled_count_100_ORF_match","assembled_count_100", "assembled_count_10_ORF_match","assembled_count_10", "not_assembled")))
ens_gene_cov$organism = "h. sapiens"
ens_gene_cov$organism[ens_gene_cov$org == "zebrafish"] = "d. rerio"
ens_gene_cov$organism[ens_gene_cov$org == "arabidopsis"] = "a. thaliana"
ens_gene_cov$organism[ens_gene_cov$org == "yeast_small"] = "s. cerevisiae"
write.table(ens_gene_cov, file="data/ensembl_pcgene_coverage.txt", sep='\t', quote=F, row.names = F)

#### BOTH STRANDS

ens_gene_cov = NULL
for(i in seq_along(orgs)){
    org = orgs[i]
    
    tt_bs = fread(paste0("data/", org, "/trinity_translated_both_strands.txt"), data.table = F)
    tt_bs$org = org
    
    load(paste0("data/", org, "/processed_ensembl.Rdata"))
    if(org == "yeast"){
        tt_bs$gene_id = ensembl_gff_cds$gene_id[match(strv_split(tt_bs$query, "[.]",1), ensembl_gff_cds$transcript_id)]
    }else if(org=="human"){
        tt_bs$gene_id = ensembl_gff_cds$gene_id[match(strv_split(tt_bs$query_gene, "[.]",1), ensembl_gff_cds$transcript_id)]
    }else if(org=="arabidopsis"){
        tt_bs$gene_id = ensembl_gff_cds$gene_id[match((tt_bs$query_gene), ensembl_gff_cds$transcript_id)]
    }else if(org=="zebrafish"){
        tt_bs$gene_id = ensembl_gff_cds$gene_id[match(strv_split(tt_bs$query_gene, "[.]",1), ensembl_gff_cds$transcript_id)]
    }
    
    ens_gene_cov_org = make_ens_coverage(org, ensembl_gff_cds, tt_bs)
    ens_gene_cov = rbind(ens_gene_cov,ens_gene_cov_org)
    
}
ens_gene_cov$class = factor(ens_gene_cov$class, levels = rev(c("assembled_count_100_ORF_match","assembled_count_100", "assembled_count_10_ORF_match","assembled_count_10", "not_assembled")))
ens_gene_cov$organism = "h. sapiens"
ens_gene_cov$organism[ens_gene_cov$org == "zebrafish"] = "d. rerio"
ens_gene_cov$organism[ens_gene_cov$org == "arabidopsis"] = "a. thaliana"
ens_gene_cov$organism[ens_gene_cov$org == "yeast_small"] = "s. cerevisiae"

write.table(ens_gene_cov, file="data/ensembl_pcgene_coverage_bothstrands.txt", sep='\t', quote=F, row.names = F)

