# combine lcs back

library(data.table)
library(tidyverse)
library(stringr)
library(halpme)
library(GenomicRanges)
library(stringdist)
library(Biostrings)

options(stringsAsFactors = F)

alt_orgs = c('human', 'human', 'yeast', 'yeast')
alt_versions = c("cds", "pc", "2sample", "3sample")

for(i in 1:length(alt_orgs)){
    org = alt_orgs[i]
    version = alt_versions[i]
    org_version = ifelse(version=="", org, paste(org, version, sep="_"))
    
    if(!file.exists(paste0("data/",org_version, "/LCS/lcs_data_summary.rdata"))){
    
        load(paste0("data/",org_version,"/ALLDATA_upto_lcs.rdata"))
        
        # first, do lcs for sequences on both strands
        # do_neg has seqs for both lcs calcs on both strands
        # do pos has above seqs AND positive strand only matches
        trinity_translate = trinity[trinity$transcript_id %in% c(e2t_best_matches_do_neg$subject),]
        e2t = e2t_best_matches_do_neg
        
        
        load(paste0("data/", org_version, "/LCS/lcs_p1.rdata"))
        load(paste0("data/", org_version, "/LCS/lcs_p2.rdata"))
        load(paste0("data/", org_version, "/LCS/lcs_p3.rdata"))  
        
        lcs_df_names = c("lcs_seq","a_start", "b_start","a_end","b_end")
        
        colnames(lcs1p) = lcs_df_names 
        colnames(lcs2p) = lcs_df_names 
        colnames(lcs3p) = lcs_df_names 
        
        z1 = as.numeric(nchar(lcs1p$lcs_seq))
        z2 = as.numeric(nchar(lcs2p$lcs_seq))
        z3 = as.numeric(nchar(lcs3p$lcs_seq))
        
        z = pmax(z1,z2,z3)
        
        zdf = data.frame(z1,z2,z3)
        zdf$lcs_frame = unlist(lapply(apply(zdf,1, function(x) which.max(x)), function(x) x[1]))
        
        zdf$lcs_len = zdf$z1
        zdf$lcs_len[which(zdf$lcs_frame == 2)] = zdf$z2[which(zdf$lcs_frame == 2)]
        zdf$lcs_len[which(zdf$lcs_frame == 3)] = zdf$z3[which(zdf$lcs_frame == 3)]
        
        
        lcs_dfp = lcs1p
        lcs_dfp[which(zdf$lcs_frame == 2),] = lcs2p[which(zdf$lcs_frame == 2),]
        lcs_dfp[which(zdf$lcs_frame == 3),] = lcs3p[which(zdf$lcs_frame == 3),]
        
        lcs_dfp$lcs_seq = as.character(unlist(lcs_dfp$lcs_seq))
        
        lcs_dfp = cbind(zdf, lcs_dfp)
        lcs_dfp$lcs_strand = "+"
        
        lcs_df = lcs_dfp
        
        lcs_df$lcs_seq = as.character(unlist(lcs_df$lcs_seq))
        for(i in 7:10){
            lcs_df[,i] = as.numeric(unlist(lcs_df[,i]))
        }
        lcs_df$first_AA = ifelse(str_sub(lcs_df$lcs_seq,1,1)=="M", "MET", "ALT")
        lcs_df$last_AA = ifelse(str_sub(lcs_df$lcs_seq,-1,-1)=="*", "STOP", "ALT")
        lcs_df$lcs_len = ifelse(lcs_df$last_AA == "STOP", lcs_df$lcs_len-1, lcs_df$lcs_len)
        
        lcs_pos_strand = lcs_df[,-c(1:3)]
        lcs_pos_strand = cbind(e2t_best_matches_do_pos[,c(1,2)], lcs_pos_strand)
        
        #lcs_pos_strand_only = lcs_pos_strand[paste0(lcs_pos_strand$subject,lcs_pos_strand$query) %in% paste0(e2t_best_matches_pos$subject, e2t_best_matches_pos$query),]
        #lcs_both_strand_only = lcs_pos_strand[paste0(lcs_pos_strand$subject,lcs_pos_strand$query) %in% paste0(e2t_best_matches$subject, e2t_best_matches$query),]
        
        load(paste0("data/", org_version, "/LCS/lcs_n1.rdata"))
        load(paste0("data/", org_version, "/LCS/lcs_n2.rdata"))
        load(paste0("data/", org_version, "/LCS/lcs_n3.rdata")) 
        
        ###########
        colnames(lcs1n) = lcs_df_names 
        colnames(lcs2n) = lcs_df_names 
        colnames(lcs3n) = lcs_df_names 
        
        z1 = as.numeric(nchar(lcs1n$lcs_seq))
        z2 = as.numeric(nchar(lcs2n$lcs_seq))
        z3 = as.numeric(nchar(lcs3n$lcs_seq))
        
        z = pmax(z1,z2,z3)
        
        zdf = data.frame(z1,z2,z3)
        zdf$lcs_frame = unlist(lapply(apply(zdf,1, function(x) which.max(x)), function(x) x[1]))
        
        zdf$lcs_len = zdf$z1
        zdf$lcs_len[which(zdf$lcs_frame == 2)] = zdf$z2[which(zdf$lcs_frame == 2)]
        zdf$lcs_len[which(zdf$lcs_frame == 3)] = zdf$z3[which(zdf$lcs_frame == 3)]
        
        
        lcs_dfn = lcs1n
        lcs_dfn[which(zdf$lcs_frame == 2),] = lcs2n[which(zdf$lcs_frame == 2),]
        lcs_dfn[which(zdf$lcs_frame == 3),] = lcs3n[which(zdf$lcs_frame == 3),]
        
        lcs_dfn$lcs_seq = as.character(unlist(lcs_dfn$lcs_seq))
        lcs_dfn = cbind(zdf, lcs_dfn)
        lcs_dfn$lcs_strand = "-"
        
        lcs_df = lcs_dfn
        
        lcs_df$lcs_seq = as.character(unlist(lcs_df$lcs_seq))
        for(i in 7:10){
            lcs_df[,i] = as.numeric(unlist(lcs_df[,i]))
        }
        lcs_df$first_AA = ifelse(str_sub(lcs_df$lcs_seq,1,1)=="M", "MET", "ALT")
        lcs_df$last_AA = ifelse(str_sub(lcs_df$lcs_seq,-1,-1)=="*", "STOP", "ALT")
        lcs_df$lcs_len = ifelse(lcs_df$last_AA == "STOP", lcs_df$lcs_len-1, lcs_df$lcs_len)
        
        lcs_neg_strand = lcs_df[,-c(1:3)]
        lcs_neg_strand = cbind(e2t_best_matches_do_neg[,c(1,2)], lcs_neg_strand)
        
        lcs_pos_strand = lcs_pos_strand[match(paste0(lcs_neg_strand$query, lcs_neg_strand$subject), 
                                              paste0(lcs_pos_strand$query, lcs_pos_strand$subject)),]
        
        pn = mapply(function(x,y) which.max(c(x,y))[1], x=lcs_pos_strand$lcs_len, y=lcs_neg_strand$lcs_len)
        lcs_df = lcs_pos_strand
        lcs_df[which(pn == 2),] = lcs_neg_strand[which(pn == 2),]
        
        lcs_both_strands = lcs_df
        
        # next, do lcs for sequences on positive strands ONLY
        # do_neg has seqs for both lcs calcs on both strands
        # do pos has above seqs AND positive strand only matches
        trinity_translate = trinity[trinity$transcript_id %in% c(e2t_best_matches_do_pos$subject),]
        e2t = e2t_best_matches_do_pos
        
        
        load(paste0("data/", org_version, "/LCS/lcs_p1.rdata"))
        load(paste0("data/", org_version, "/LCS/lcs_p2.rdata"))
        load(paste0("data/", org_version, "/LCS/lcs_p3.rdata"))  
        
        lcs_df_names = c("lcs_seq","a_start", "b_start","a_end","b_end")
        
        colnames(lcs1p) = lcs_df_names 
        colnames(lcs2p) = lcs_df_names 
        colnames(lcs3p) = lcs_df_names 
        
        z1 = as.numeric(nchar(lcs1p$lcs_seq))
        z2 = as.numeric(nchar(lcs2p$lcs_seq))
        z3 = as.numeric(nchar(lcs3p$lcs_seq))
        
        z = pmax(z1,z2,z3)
        
        zdf = data.frame(z1,z2,z3)
        zdf$lcs_frame = unlist(lapply(apply(zdf,1, function(x) which.max(x)), function(x) x[1]))
        
        zdf$lcs_len = zdf$z1
        zdf$lcs_len[which(zdf$lcs_frame == 2)] = zdf$z2[which(zdf$lcs_frame == 2)]
        zdf$lcs_len[which(zdf$lcs_frame == 3)] = zdf$z3[which(zdf$lcs_frame == 3)]
        
        
        lcs_dfp = lcs1p
        lcs_dfp[which(zdf$lcs_frame == 2),] = lcs2p[which(zdf$lcs_frame == 2),]
        lcs_dfp[which(zdf$lcs_frame == 3),] = lcs3p[which(zdf$lcs_frame == 3),]
        
        lcs_dfp$lcs_seq = as.character(unlist(lcs_dfp$lcs_seq))
        
        lcs_dfp = cbind(zdf, lcs_dfp)
        lcs_dfp$lcs_strand = "+"
        
        lcs_df = lcs_dfp
        
        lcs_df$lcs_seq = as.character(unlist(lcs_df$lcs_seq))
        for(i in 7:10){
            lcs_df[,i] = as.numeric(unlist(lcs_df[,i]))
        }
        lcs_df$first_AA = ifelse(str_sub(lcs_df$lcs_seq,1,1)=="M", "MET", "ALT")
        lcs_df$last_AA = ifelse(str_sub(lcs_df$lcs_seq,-1,-1)=="*", "STOP", "ALT")
        lcs_df$lcs_len = ifelse(lcs_df$last_AA == "STOP", lcs_df$lcs_len-1, lcs_df$lcs_len)
        
        lcs_pos_strand = lcs_df[,-c(1:3)]
        lcs_pos_strand = cbind(e2t_best_matches_do_pos[,c(1,2)], lcs_pos_strand)
        
        # still need to filter out 'pos only' as this contains 'both strands' as well...
        
        save(lcs_both_strands, lcs_pos_strand, file=paste0("data/",org_version, "/LCS/lcs_data_summary.rdata"))
    }
}
