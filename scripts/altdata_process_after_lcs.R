## BOTH STRANDS

library(data.table)
library(tidyverse)
library(stringr)
library(halpme)
library(GenomicRanges)
library(stringdist)
library(Biostrings)

options(stringsAsFactors = F)
source("scripts/helper_functions.R")
org="yeast"



alt_orgs = c('human', 'human', 'yeast', 'yeast')
alt_versions = c("cds", "pc", "2sample", "3sample")

for(i in 1:length(alt_orgs)){
    org = alt_orgs[i]
    version = alt_versions[i]
    org_version = ifelse(version=="", org, paste(org, version, sep="_"))
    
    
    if(!file.exists(paste0('data/', org_version, "/trinity_translated_both_strands.txt")) | !file.exists(paste0('data/', org_version, "/trinity_translated.txt"))){
    
        for(do_set in c("both_strands", "pos_strand")){
        
            load(paste0("data/", org_version, "/ALLDATA_upto_lcs.rdata"))
            load(paste0("data/", org, "/processed_ensembl.Rdata"))
            
            ## commented code from (process_upto_lcs.R)
            # e2t_best_matches_do_pos = rbind(e2t_best_matches, e2t_best_matches_pos) %>% distinct()
            # e2t_best_matches_do_neg = e2t_best_matches
            # 
            # save(e2t_best_matches_do_pos, e2t_best_matches_do_neg, trinity, ensembl_peps,longest_common_substring, file=paste0("data/", org, "/upto_lcs.rdata"))
            # save.image(paste0("data/",org, "/ALLDATA_upto_lcs.rdata"))
            
            # e2t_best_matches has pos + neg
            # e2t_best_matches_pos has positive strand only
            
            if(do_set == "both_strands"){
                trinity_translate = trinity[trinity$transcript_id %in% c(e2t_best_matches$subject),]
            }else{
                trinity_translate = trinity[trinity$transcript_id %in% c(e2t_best_matches_do_pos$subject),]
            }
            #lcs_both_strands_v2 = cbind(e2t_best_matches_do_neg_v2[,c(1,2)], lcs_both)
            
            load(paste0("data/",org_version, "/LCS/lcs_data_summary.rdata"))
            
            # LCS for both strands
            if(do_set == "both_strands"){
                longest_common_subseq_borf = lcs_both_strands
            }else{
                longest_common_subseq_borf = lcs_pos_strand
            }
            colnames(longest_common_subseq_borf)[c(6,7,8,9,11,12)] = c("lcs_trinity_start","lcs_ensembl_start", "lcs_trinity_end","lcs_ensembl_end","lcs_first_AA", "lcs_last_AA")
            
            longest_common_subseq_borf$lcs_start_site_nt = (longest_common_subseq_borf$lcs_trinity_start*3) -3 + longest_common_subseq_borf$lcs_frame
            longest_common_subseq_borf$lcs_stop_site_nt = (longest_common_subseq_borf$lcs_len*3) + longest_common_subseq_borf$lcs_start_site_nt -1
            longest_common_subseq_borf$lcs_utr3_length = (trinity$seq_len[match(longest_common_subseq_borf$subject, trinity$transcript_id)] - 
                                                              longest_common_subseq_borf$lcs_stop_site_nt)
            longest_common_subseq_borf$trinity_length_nt = trinity$seq_len[match(longest_common_subseq_borf$subject, trinity$transcript_id)]
            
            #### add in trinity transcript classes ####
            # NOTE: partial:: usually splicing variants
            
            if(do_set == "both_strands"){
                single_match = left_join(e2t_single_gene_match, longest_common_subseq_borf, by=c('query'='query', 'subject'='subject'), suffix = c('', ".lcs"))
                multi_match = left_join(e2t_multi_gene_match, longest_common_subseq_borf, by=c('query'='query', 'subject'='subject'), suffix = c('', ".lcs"))
            }else{
                single_match = left_join(e2t_single_gene_match_pos, longest_common_subseq_borf, by=c('query'='query', 'subject'='subject'), suffix = c('', ".lcs"))
                multi_match = left_join(e2t_multi_gene_match_pos, longest_common_subseq_borf, by=c('query'='query', 'subject'='subject'), suffix = c('', ".lcs"))
            }
            single_match$match_type = 'single_hit'
            single_match = single_match[paste0(single_match$query, single_match$subject) %in% paste0(longest_common_subseq_borf$query, longest_common_subseq_borf$subject),]
            multi_match$match_type = 'multi_hit'
            multi_match = multi_match[paste0(multi_match$query, multi_match$subject) %in% paste0(longest_common_subseq_borf$query, longest_common_subseq_borf$subject),]
            trinity_translated = plyr::rbind.fill(single_match, multi_match)
            trinity_translated$ens_aa_len = nchar(ensembl_peps$seq[match(trinity_translated$query, ensembl_peps$transcript)])
            
            trinity_translated$base_class = "insufficient_coverage"
            
            trinity_translated$base_class[trinity_translated$lcs_len >=10 & trinity_translated$lcs_first_AA == "ALT" & trinity_translated$lcs_last_AA == "ALT" ] = "incomplete"
            trinity_translated$base_class[trinity_translated$lcs_len >=10 & trinity_translated$lcs_first_AA == "ALT" & trinity_translated$lcs_last_AA == "ALT" & 
                                              (trinity_translated$lcs_trinity_end - (trinity_translated$trinity_length_nt/3) < -2 | trinity_translated$lcs_trinity_start > 1)] = "incomplete_partial"
            
            trinity_translated$base_class[trinity_translated$lcs_len >=10 & trinity_translated$lcs_first_AA == "MET" & trinity_translated$lcs_last_AA == "ALT" & trinity_translated$lcs_trinity_end - (trinity_translated$trinity_length_nt/3) >= -2] = "incomplete_3prime"
            trinity_translated$base_class[trinity_translated$lcs_len >=10 & trinity_translated$lcs_first_AA == "MET" & trinity_translated$lcs_last_AA == "ALT" & trinity_translated$lcs_trinity_end - (trinity_translated$trinity_length_nt/3) < -2] = "incomplete_3prime_partial"
            
            trinity_translated$base_class[trinity_translated$lcs_len >=10 & trinity_translated$lcs_first_AA == "ALT" & trinity_translated$lcs_last_AA == "STOP" & trinity_translated$lcs_trinity_start == 1] = "incomplete_5prime"
            trinity_translated$base_class[trinity_translated$lcs_len >=10 & trinity_translated$lcs_first_AA == "ALT" & trinity_translated$lcs_last_AA == "STOP" & trinity_translated$lcs_trinity_start > 1] = "incomplete_5prime_partial"
            
            
            trinity_translated$base_class[trinity_translated$lcs_len >=10 & trinity_translated$lcs_first_AA == "MET" & trinity_translated$lcs_last_AA == "STOP" ] = "complete_partial"
            trinity_translated$base_class[trinity_translated$lcs_len >=10 & trinity_translated$lcs_first_AA == "MET" & trinity_translated$lcs_last_AA == "STOP" & trinity_translated$lcs_len == trinity_translated$ens_aa_len] = "complete"
            
            ####### get ensembl M location relative to trinity start (in NT)
            
            # if base_class == complete/ic 3prime
            #incomplete/ic 5 prime - b_start*3+frame-1
            
            trinity_translated$dist_to_ens_M = NA
            trinity_translated$dist_to_ens_M[trinity_translated$base_class %in% c("complete", "complete_partial", "incomplete_3prime","incomplete_3prime_partial")] = trinity_translated$lcs_start_site_nt[trinity_translated$base_class %in% c("complete", "complete_partial", "incomplete_3prime","incomplete_3prime_partial")]
            index = which(trinity_translated$base_class %in% c("incomplete", "incomplete_partial", "incomplete_5prime","incomplete_5prime_partial"))
            trinity_translated$dist_to_ens_M[index] = (trinity_translated$lcs_ensembl_start[index]*3+trinity_translated$lcs_frame[index]-1)*-1
            
            trinity_seqs = trinity[match(trinity_translated$subject, trinity$transcript_id),]
            trinity_seqs$frame = trinity_translated$lcs_frame
            trinity_seqs$translated = NA
            
            get_translated = function(trinity_seqs, frame=1, strand="+"){
                
                if(strand == "+"){
                trinity_seqs$translated[which(trinity_seqs$frame == frame & trinity_translated$lcs_strand == strand)] = 
                    as.character(translate(DNAStringSet(str_sub(trinity_seqs$seq[which(trinity_seqs$frame == frame & trinity_translated$lcs_strand == strand)],1,-1))))
                }else{
                    trinity_seqs$translated[which(trinity_seqs$frame == frame & trinity_translated$lcs_strand == strand)] = 
                        as.character(translate(reverseComplement(DNAStringSet(str_sub(trinity_seqs$seq[which(trinity_seqs$frame == 1 & trinity_translated$lcs_strand == strand)],1,-1*(frame))))))
                }
                
                return(trinity_seqs)
            }
            
            if(do_set == "both_strands"){
                for(strand in c("+", "-")){
                    for(frame in c(1:3)){
                        trinity_seqs = get_translated(trinity_seqs=trinity_seqs, frame=frame, strand=strand)
                    }
                }
            }else{
                for(frame in c(1:3)){
                    trinity_seqs = get_translated(trinity_seqs=trinity_seqs, frame=frame, strand="+")
                }
            }
            
            trinity_seqs$start_loc = trinity_translated$lcs_trinity_start
            upstream_seqs = str_sub(trinity_seqs$translated, 1, trinity_seqs$start_loc-1)
            upstream_stop = str_locate(stringi::stri_reverse(upstream_seqs), "[*]")[,1]
            upstream_stop[which(!grepl("[*]", upstream_seqs))] = nchar(upstream_seqs[which(!grepl("[*]", upstream_seqs))])
            
            summary(upstream_stop)
            trinity_translated$upstream_stop = upstream_stop
            trinity_translated$has_upstream_stop = ifelse(grepl("[*]", upstream_seqs), 1,0)
            
            trinity_translated$first_M = str_locate(trinity_seqs$translated, "M")[,1]
            aggregate(first_M ~ base_class, trinity_translated[trinity_translated$lcs_strand == "+" & trinity_translated$base_class %in% c("incomplete", "incomplete_5prime"),], summary)
            
            table(trinity_translated$base_class)
            #trinity_translated$seq_sim = stringsim(trinity_seqs$translated, ensembl_peps$seq[match(trinity_translated$query, ensembl_peps$transcript)])
            trinity_translated$lcs_seq_sim = stringsim(gsub("[*]","",trinity_translated$lcs_seq), ensembl_peps$seq[match(trinity_translated$query, ensembl_peps$transcript)])
            trinity_translated$lcs_percent_of_ens = 1- (trinity_translated$ens_aa_len - trinity_translated$lcs_len) / trinity_translated$ens_aa_len
            trinity_translated$lcs_diff_to_ens =(trinity_translated$ens_aa_len - trinity_translated$lcs_len) 
            
            #### UTR coverage ####
            
            # 5 'utr_len = dist_to_ensM 
            # average 293nt (in tx with 5' end of ORF)
            summary(trinity_translated$dist_to_ens_M[trinity_translated$base_class %in% c("complete", "incomplete_3prime")])
            summary(trinity_translated$lcs_start_site_nt[trinity_translated$base_class %in% c("complete", "incomplete_3prime")])
            # 3' UTR len average 281nt
            summary(trinity_translated$lcs_utr3_length[trinity_translated$base_class %in% c("complete", "incomplete_5prime")])
            utr3_index = which(tolower(utr_annot$type)=="three_prime_utr")
            if(org == "human" | org == "zebrafish"){ 
                tt_index = match(strv_split(trinity_translated$query, "[.]",1), utr_annot$transcript[utr3_index])
            }else{
                tt_index = match(trinity_translated$query, utr_annot$transcript[utr3_index])
            }
            trinity_translated$ens_3utr_len = utr_annot$end[utr3_index][tt_index] - utr_annot$start[utr3_index][tt_index]
            utr5_index = which(tolower(utr_annot$type)=="five_prime_utr")
            if(org == "human" | org == "zebrafish"){
                tt_index = match(strv_split(trinity_translated$query, "[.]",1), utr_annot$transcript[utr5_index])
            }else{
                tt_index = match(trinity_translated$query, utr_annot$transcript[utr5_index])
            }
            trinity_translated$ens_5utr_len = utr_annot$end[utr5_index][tt_index] - utr_annot$start[utr5_index][tt_index]
            
            summary(trinity_translated$ens_3utr_len)
            summary(trinity_translated$ens_5utr_len)
            if(do_set == "both_strands"){
                write.table(trinity_translated, paste0('data/', org_version, "/trinity_translated_both_strand_partials.txt"), sep='\t', quote=F, row.names = F)
            }else{
                write.table(trinity_translated, paste0('data/', org_version, "/trinity_translated_partials.txt"), sep='\t', quote=F, row.names = F)
            }
            # extend 'partial' matches and recalculate some stuff...
            
            # add columns for TRINITY ORFs (not just LCS)
            trinity_translated$trinity_orf_start = trinity_translated$lcs_trinity_start
            trinity_translated$trinity_orf_end = trinity_translated$lcs_trinity_end
            trinity_translated$trinity_orf_length = trinity_translated$lcs_trinity_end - trinity_translated$lcs_trinity_start
            trinity_translated$trinity_start_site_nt = (trinity_translated$trinity_orf_start*3) -3 + trinity_translated$lcs_frame
            trinity_translated$trinity_stop_site_nt = (trinity_translated$trinity_orf_length*3) + trinity_translated$trinity_start_site_nt -1
            trinity_translated$trinity_utr3_length = (trinity$seq_len[match(trinity_translated$subject, trinity$transcript_id)] - 
                                                          trinity_translated$trinity_stop_site_nt)
            trinity_translated$trinity_length_nt = trinity$seq_len[match(trinity_translated$subject, trinity$transcript_id)]
            trinity_translated$trinity_first_AA = trinity_translated$lcs_first_AA
            trinity_translated$trinity_last_AA = trinity_translated$lcs_last_AA
            trinity_translated$trinity_seqsim_ensembl = NA
            
            ## INCOMPLETE 3' PARTIAL MATCHES
            
            ## INCOMPLETE 3' PARTIAL MATCHES
            tt_3pp = trinity_translated[which(trinity_translated$base_class == "incomplete_3prime_partial"),]
            trinity_translate = trinity[match(tt_3pp$subject, trinity$transcript_id),]
            
            translate_by_frame = function(seqs, frames, strands){
                
                aa1p = translate(DNAStringSet(seqs))
                aa2p = translate(DNAStringSet(str_sub(seqs,2,-1)))
                aa3p = translate(DNAStringSet(str_sub(seqs,3,-1)))
                aa1n = translate(reverseComplement(DNAStringSet(seqs)))
                aa2n = translate(reverseComplement(DNAStringSet(str_sub(seqs,1,-2))))
                aa3n = translate(reverseComplement(DNAStringSet(str_sub(seqs,1,-3))))
                
                aa = aa1p
                aa[frames == 2 & strands == "+"] = aa2p[frames == 2 & strands == "+"]
                aa[frames == 3 & strands == "+"] = aa3p[frames == 3 & strands == "+"]
                aa[frames == 1 & strands == "-"] = aa1n[frames == 1 & strands == "-"]
                aa[frames == 2 & strands == "-"] = aa2n[frames == 2 & strands == "-"]
                aa[frames == 3 & strands == "-"] = aa3n[frames == 3 & strands == "-"]
                
                aa = as.character(aa)
                return(aa)
                
            }
            
            translate_by_frame_pos = function(seqs, frames){
                
                aa1p = translate(DNAStringSet(seqs))
                aa2p = translate(DNAStringSet(str_sub(seqs,2,-1)))
                aa3p = translate(DNAStringSet(str_sub(seqs,3,-1)))
                   
                aa = aa1p
                aa[frames == 2] = aa2p[frames == 2]
                aa[frames == 3] = aa3p[frames == 3]
                
                aa = as.character(aa)
                return(aa)
                
            }
            if(do_set == "both_strands"){
                aa = translate_by_frame(trinity_translate$seq, frames = tt_3pp$lcs_frame, strands = tt_3pp$match_strand) 
            }else{
                aa = translate_by_frame_pos(trinity_translate$seq, frames = tt_3pp$lcs_frame) 
            }
            # FULL AA sequence -- from the start of the LCS (M, as 5' should be complete), to the end
            full_aa = mapply(function(x,y) str_sub(x, (str_locate(x, y)[1]), -1), x=aa, y=tt_3pp$lcs_seq)
            full_aa = as.character(full_aa)
            stop_locs = str_locate(full_aa, "[*]")[,2]
            
            #tt_3pp$lcs_trinity_start = tt_3pp$lcs_trinity_start
            
            # sequence upstream of the LCS (M)
            upstream_aa = as.character(mapply(function(x,y) str_sub(x, 1, (str_locate(gsub("[*]","Z",x), gsub("[*]","Z",y))[1] -1)), x=aa, y=full_aa))
            tt_3pp$upstream_aa = upstream_aa
            
            stop_locs_for_orf = stop_locs
            stop_locs_for_orf[is.na(stop_locs_for_orf)] = nchar(full_aa)[is.na(stop_locs_for_orf)]
            
            # from LCS M to the next * ... i.e. extended the incomplete 3' end
            trinity_orf = as.character(mapply(function(x,y) str_sub(x, 1, y), x=full_aa, y=stop_locs_for_orf))
            
            # get locations back again
            aa_locs = mapply(function(x,y) str_locate(x, y), x=aa, y=trinity_orf)
            dimnames(aa_locs) = NULL
            aa_locs = t(aa_locs)
            
            tt_3pp$trinity_orf_start = aa_locs[,1]
            tt_3pp$trinity_orf_end = aa_locs[,2]
            tt_3pp$trinity_last_AA[!is.na(stop_locs)] = "STOP"
            
            tt_3pp$trinity_orf_length = tt_3pp$trinity_orf_end - tt_3pp$trinity_orf_start
            tt_3pp$trinity_start_site_nt = (tt_3pp$trinity_orf_start*3) -3 + tt_3pp$lcs_frame
            tt_3pp$trinity_stop_site_nt = (tt_3pp$trinity_orf_length*3) + tt_3pp$trinity_start_site_nt -1
            tt_3pp$trinity_utr3_length = (trinity$seq_len[match(tt_3pp$subject, trinity$transcript_id)] - 
                                              tt_3pp$trinity_stop_site_nt)
            
            
            # where is the new STOP?
            aa_rev = gsub("[*]","Z", stringi::stri_reverse(aa))
            orf_seq_rev = gsub("[*]","Z", stringi::stri_reverse(trinity_orf))
            end_locs = as.numeric(mapply(function(x,y) str_locate(x, y)[1], x=aa_rev, y=orf_seq_rev))
            
            
            # redo the base classes
            tt_3pp$base_class = "incomplete_3prime_partial"
            # no stop in extended seq... still incomplete 3'
            tt_3pp$base_class[end_locs == 1 & tt_3pp$trinity_orf_length >= tt_3pp$lcs_len] = "incomplete_3prime"
            # stop in extended seq & longer than the initial LCS
            tt_3pp$base_class[tt_3pp$trinity_last_AA == "STOP" & tt_3pp$trinity_orf_end >= tt_3pp$lcs_trinity_end & 
                                  tt_3pp$trinity_orf_start <= tt_3pp$lcs_trinity_start] = "complete"
            
            # seq_similiarity to ensembl pep
            tt_3pp$trinity_seqsim_ensembl = stringsim(gsub("[*]","",trinity_orf), 
                                                      ensembl_peps$seq[match(tt_3pp$query, ensembl_peps$transcript)])
            tt_3pp$trinity_orf = trinity_orf
            
            
            #######################################
            
            
            ## INCOMPLETE PARTIAL 5' MATCHES
            
            tt_5pp = trinity_translated[trinity_translated$base_class == "incomplete_5prime_partial",]
            trinity_translate = trinity[match(tt_5pp$subject, trinity$transcript_id),]
            
            if(do_set == "both_strands"){
                aa = translate_by_frame(trinity_translate$seq, frames = tt_5pp$lcs_frame, strands = tt_5pp$match_strand) 
            }else{
                aa = translate_by_frame_pos(trinity_translate$seq, frames = tt_5pp$lcs_frame) 
            }
            
            # reverse & find first stop after the start 'M' / start location (5' IC)
            aa_rev = gsub("[*]","Z", stringi::stri_reverse(aa))
            lcs_seq_rev = gsub("[*]","Z", stringi::stri_reverse(tt_5pp$lcs_seq))
            full_aa = mapply(function(x,y) str_sub(x, (str_locate(x, y)[1]), -1), x=aa_rev, y=lcs_seq_rev)
            full_aa = str_sub(as.character(full_aa),2,-1)
            stop_locs = str_locate(full_aa, "Z")[,2]
            stop_locs[is.na(stop_locs)] = nchar(full_aa)[is.na(stop_locs)]
            full_aa = as.character(mapply(function(x,y) str_sub(x, 1, y-1),x=full_aa, y=stop_locs))
            # reverse again (back to normal direction)
            full_aa = stringi::stri_reverse(full_aa)
            # and find any 'M's
            start_locs = str_locate(full_aa, "M")[,2]
            start_locs[is.na(start_locs)] = 1
            # complete ORF sequence
            full_aa = as.character(mapply(function(x,y) str_sub(x, y, -1),x=full_aa, y=start_locs))
            aa_locs = mapply(function(x,y) str_locate(x, y), x=aa, y=full_aa)
            dimnames(aa_locs) = NULL
            aa_locs = t(aa_locs)
            
            tt_5pp$trinity_orf_start = aa_locs[,1]
            tt_5pp$trinity_orf_end = aa_locs[,2] +1
            tt_5pp$trinity_first_AA[str_sub(full_aa, 1,1) == "M"] ="MET"
            tt_5pp$base_class = "incomplete_5prime"
            tt_5pp$base_class[tt_5pp$trinity_first_AA == "MET"] = "complete"
            # is the LCS start upstream of the ORF start?
            is_shorter = tt_5pp$lcs_trinity_start < aa_locs[,1]
            tt_5pp$base_class[tt_5pp$trinity_first_AA == "MET" & is_shorter] = "complete_partial" # changed from 'complete partial'
            upstream_aa = as.character(mapply(function(x,y) str_sub(x, 1, (str_locate(gsub("[*]","Z",x), gsub("[*]","Z",y))[1] -1)), x=aa, y=full_aa))
            tt_5pp$upstream_aa = upstream_aa
            
            tt_5pp$trinity_orf_length = tt_5pp$trinity_orf_end - tt_5pp$trinity_orf_start
            tt_5pp$trinity_start_site_nt = (tt_5pp$trinity_orf_start*3) -3 + tt_5pp$lcs_frame
            tt_5pp$trinity_stop_site_nt = (tt_5pp$trinity_orf_length*3) + tt_5pp$trinity_start_site_nt -1
            tt_5pp$trinity_utr3_length = (trinity$seq_len[match(tt_5pp$subject, trinity$transcript_id)] - 
                                              tt_5pp$trinity_stop_site_nt)
            trinity_orf = full_aa
            # seq_similiarity to ensembl pep?
            tt_5pp$trinity_seqsim_ensembl = stringsim(gsub("[*]","",trinity_orf), 
                                                      ensembl_peps$seq[match(tt_5pp$query, ensembl_peps$transcript)])
            tt_5pp$trinity_orf = trinity_orf
            
            
            ## INCOMPLETE MATCHES (3' AND 5')
            
            tt_ic = trinity_translated[trinity_translated$base_class == "incomplete_partial",]
            trinity_translate = trinity[match(tt_ic$subject, trinity$transcript_id),]
            
            if(do_set == "both_strands"){
                aa = translate_by_frame(trinity_translate$seq, frames = tt_ic$lcs_frame, strands = tt_ic$match_strand) 
            }else{
                aa = translate_by_frame_pos(trinity_translate$seq, frames = tt_ic$lcs_frame) 
            }
            
            # check for anything past the end (3' ic-ness)
            full_aa = mapply(function(x,y) str_sub(x, (str_locate(gsub("[*]","Z",x), gsub("[*]","Z",y))[1]), -1), x=aa, y=tt_ic$lcs_seq)
            full_aa = as.character(full_aa)
            #stop_locs = str_locate(full_aa, "[*]")[,2]
            # check for anythin before start (5' ic-ness)
            aa_rev = gsub("[*]","Z", stringi::stri_reverse(aa))
            lcs_seq_rev = gsub("[*]","Z", stringi::stri_reverse(tt_ic$lcs_seq))
            full_aa = mapply(function(x,y) str_sub(x, (str_locate(x, y)[1]), -1), x=aa_rev, y=lcs_seq_rev)
            full_aa = as.character(str_sub(as.character(full_aa),2,-1))
            stop_locs = str_locate(full_aa, "Z")[,2]
            stop_locs[is.na(stop_locs)] = nchar(full_aa)[is.na(stop_locs)]
            full_aa = as.character(mapply(function(x,y) str_sub(x, 1, y),x=full_aa, y=stop_locs))
            
            full_aa = stringi::stri_reverse(full_aa)
            start_locs = str_locate(full_aa, "M")[,2]
            start_locs[is.na(start_locs)] = 1
            full_aa = as.character(mapply(function(x,y) str_sub(x, y, -1),x=full_aa, y=start_locs))
            full_aa = mapply(function(x,y) str_sub(x, (str_locate(x, y)[1]), -1), x=aa, y=full_aa)
            full_aa = as.character(full_aa)
            # check for anything past the end (3' IC)
            stop_locs = str_locate(full_aa, "[*]")[,2]
            stop_locs[is.na(stop_locs)] = nchar(full_aa)[is.na(stop_locs)]
            full_aa = mapply(function(x,y) str_sub(x, 1, y), x=full_aa, y=stop_locs)
            full_aa = as.character(full_aa)
            full_ic_locs = mapply(function(x,y) str_locate(string = x, pattern = y), x=aa, y=full_aa)
            dimnames(full_ic_locs) = NULL
            full_ic_locs = t(full_ic_locs)
            
            
            tt_ic$trinity_orf_start = full_ic_locs[,1]
            tt_ic$trinity_orf_end = full_ic_locs[,2]
            tt_ic$trinity_last_AA[str_sub(full_aa,-1,-1) == "*"] = "STOP"
            tt_ic$trinity_first_AA[str_sub(full_aa,1,1) == "M"] = "MET"
            
            tt_ic$trinity_orf_length = tt_ic$trinity_orf_end - tt_ic$trinity_orf_start
            tt_ic$trinity_start_site_nt = (tt_ic$trinity_orf_start*3) -3 + tt_ic$lcs_frame
            tt_ic$trinity_stop_site_nt = (tt_ic$trinity_orf_length*3) + tt_ic$trinity_start_site_nt -1
            tt_ic$trinity_utr3_length = (trinity$seq_len[match(tt_ic$subject, trinity$transcript_id)] - 
                                             tt_ic$trinity_stop_site_nt)
            trinity_orf = full_aa
            
            aa_rev = gsub("[*]","Z", stringi::stri_reverse(aa))
            orf_seq_rev = gsub("[*]","Z", stringi::stri_reverse(trinity_orf))
            end_locs = as.numeric(mapply(function(x,y) str_locate(x, y)[1], x=aa_rev, y=orf_seq_rev))
            
            # seq_similiarity to ensembl pep?
            tt_ic$trinity_seqsim_ensembl = stringsim(gsub("[*]","",trinity_orf), 
                                                     ensembl_peps$seq[match(tt_ic$query, ensembl_peps$transcript)])
            
            tt_ic$base_class[tt_ic$lcs_trinity_start >= full_ic_locs[,1] & 
                                 tt_ic$lcs_trinity_end <= full_ic_locs[,2] & 
                                 str_sub(full_aa,1,1) == "M" & str_sub(full_aa,-1,-1) == "*"] = "complete" # complete, but longer than match
            tt_ic$base_class[tt_ic$lcs_trinity_start >= full_ic_locs[,1] & 
                                 tt_ic$lcs_trinity_end <= full_ic_locs[,2] & 
                                 str_sub(full_aa,1,1) == "M" & str_sub(full_aa,-1,-1) != "*"] = "incomplete_3prime_partial" # incomplete, longer than match, last AA may be wrong
            tt_ic$base_class[tt_ic$lcs_trinity_start >= full_ic_locs[,1] & 
                                 tt_ic$lcs_trinity_end <= full_ic_locs[,2] & 
                                 str_sub(full_aa,1,1) != "M" & str_sub(full_aa,-1,-1) == "*"] = "incomplete_5prime_partial" # incomplete, longer than match, first AA may be wrong
            tt_ic$base_class[end_locs == 1 & tt_ic$lcs_trinity_start >= full_ic_locs[,1] & 
                                 tt_ic$lcs_trinity_end <= full_ic_locs[,2] & str_sub(full_aa,1,1) == "M" &
                                 str_sub(full_aa,-1,-1) != "*"] = "incomplete_3prime"
            tt_ic$base_class[tt_ic$trinity_orf_start == 1 & tt_ic$lcs_trinity_start >= full_ic_locs[,1] & 
                                 tt_ic$lcs_trinity_end <= full_ic_locs[,2] & str_sub(full_aa,1,1) != "M" & 
                                 str_sub(full_aa,-1,-1) == "*"] = "incomplete_5prime"
            
            upstream_aa = as.character(mapply(function(x,y) str_sub(x, 1, (str_locate(gsub("[*]","Z",x), gsub("[*]","Z",y))[1] -1)), x=aa, y=full_aa))
            tt_ic$upstream_aa = upstream_aa
            
            tt_ic$trinity_orf = trinity_orf
            
            
            trinity_translated_np = trinity_translated[trinity_translated$base_class %in% c("complete", "incomplete",'incomplete_5prime','incomplete_3prime', 'insufficient_coverage','complete_partial'),]
            trinity_translated_np$base_class_tag = "full_match"
            trinity_translated_np$base_class_tag[trinity_translated_np$base_class == 'complete_partial'] = "partial_match"
            #trinity_translated_np$base_class[trinity_translated_np$base_class == "complete_partial"] = 'complete'
            
            tt_ic$base_class_tag = "partial_match"
            tt_5pp$base_class_tag = "partial_match"
            tt_3pp$base_class_tag = "partial_match"
            tt_part = rbind(tt_ic,tt_3pp,tt_5pp)
            upstream_stop = str_locate(stringi::stri_reverse(tt_part$upstream_aa), "[*]")[,1]
            upstream_stop[which(!grepl("[*]", tt_part$upstream_aa))] = nchar(tt_part$upstream_aa[which(!grepl("[*]", tt_part$upstream_aa))])
            tt_part$trinity_upstream_stop = upstream_stop
            tt_part$trinity_has_upstream_stop = ifelse(grepl("[*]", tt_part$upstream_aa), 1,0)
            tt_part$upstream_aa = NULL
            
            trinity_translated_np$trinity_orf = trinity_translated_np$lcs_seq
            trinity_translated_np$trinity_upstream_stop = trinity_translated_np$upstream_stop
            trinity_translated_np$trinity_has_upstream_stop = trinity_translated_np$has_upstream_stop
            
            trinity_translated_np = rbind(trinity_translated_np[,match(colnames(tt_part), colnames(trinity_translated_np))], tt_part)
            
            table(trinity_translated_np$base_class, trinity_translated_np$base_class_tag)  
            
            
            
            #table(tt_part$base_class, tt_part$base_class_tag)
            if(do_set == "both_strands"){
                write.table(trinity_translated_np, paste0('data/', org_version, "/trinity_translated_both_strands.txt"), quote=F, row.names = F, sep='\t')
            }else{
                write.table(trinity_translated_np, paste0('data/', org_version, "/trinity_translated.txt"), quote=F, row.names = F, sep='\t')
            }
        }
        if(org == "human"){
            
            load(paste0("data/", org, "/processed_ensembl.Rdata"))
            
            tt = fread(paste0('data/', org_version, "/trinity_translated.txt"), data.table = F)
            utr_annot$transcript = gsub(";","", utr_annot$transcript)
            utr3_index = which(tolower(utr_annot$type)=="three_prime_utr")
            if(org == "human" | org == "zebrafish"){ 
                tt_index = match(strv_split(tt$query, "[.]",1), utr_annot$transcript[utr3_index])
            }else{
                tt_index = match(tt$query, utr_annot$transcript[utr3_index])
            }
            tt$ens_3utr_len = utr_annot$end[utr3_index][tt_index] - utr_annot$start[utr3_index][tt_index]
            utr5_index = which(tolower(utr_annot$type)=="five_prime_utr")
            if(org == "human" | org == "zebrafish"){
                tt_index = match(strv_split(tt$query, "[.]",1), utr_annot$transcript[utr5_index])
            }else{
                tt_index = match(tt$query, utr_annot$transcript[utr5_index])
            }
            tt$ens_5utr_len = utr_annot$end[utr5_index][tt_index] - utr_annot$start[utr5_index][tt_index]
            write.table(tt, paste0('data/', org_version, "/trinity_translated.txt"), quote=F, row.names = F, sep='\t')
            
            tt = fread(paste0('data/', org_version, "/trinity_translated_both_strands.txt"), data.table = F)
            utr_annot$transcript = gsub(";","", utr_annot$transcript)
            utr3_index = which(tolower(utr_annot$type)=="three_prime_utr")
            if(org == "human" | org == "zebrafish"){ 
                tt_index = match(strv_split(tt$query, "[.]",1), utr_annot$transcript[utr3_index])
            }else{
                tt_index = match(tt$query, utr_annot$transcript[utr3_index])
            }
            tt$ens_3utr_len = utr_annot$end[utr3_index][tt_index] - utr_annot$start[utr3_index][tt_index]
            utr5_index = which(tolower(utr_annot$type)=="five_prime_utr")
            if(org == "human" | org == "zebrafish"){
                tt_index = match(strv_split(tt$query, "[.]",1), utr_annot$transcript[utr5_index])
            }else{
                tt_index = match(tt$query, utr_annot$transcript[utr5_index])
            }
            tt$ens_5utr_len = utr_annot$end[utr5_index][tt_index] - utr_annot$start[utr5_index][tt_index]
            write.table(tt, paste0('data/', org_version, "/trinity_translated_both_strands.txt"), quote=F, row.names = F, sep='\t')
            
        }
    }
    
}
