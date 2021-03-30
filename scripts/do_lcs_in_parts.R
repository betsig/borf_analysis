# lcs processing

library(data.table)
library(tidyverse)
library(stringr)
library(halpme)
library(GenomicRanges)
library(stringdist)
library(Biostrings)
library(progress)

options(stringsAsFactors = F)


args = commandArgs(trailingOnly=TRUE)
do_lcs_frame = args[1]

load(paste0("upto_lcs.rdata"))


if(str_sub(do_lcs_frame, 1,1) == "p"){
    trinity_translate = trinity[trinity$transcript_id %in% e2t_best_matches_do_pos$subject,]
    e2t = e2t_best_matches_do_pos #both pos strand only AND negative strand data
}else{
    trinity_translate = trinity[trinity$transcript_id %in% e2t_best_matches_do_neg$subject,]
    e2t = e2t_best_matches_do_neg
}

chunk_size = 10
strand = "both"

make_lcs_df_chunk = function(chunk_size, trinity_aa, ensembl_aa){
    
    chunks = ceiling(length(trinity_aa) / chunk_size)
    lcs = NULL
    pbj = progress_bar$new(total = chunks, format = "[:bar] :current/:total (:percent) time: :elapsedfull eta: :eta")
    for(i in 1:chunks){
        
        start = i*chunk_size - chunk_size + 1
        end = min(length(trinity_aa), start+chunk_size-1)
        lcs_result = 
            
            lcs = rbind(lcs, matrix(mapply(function(x,y) longest_common_substring(x, y, return="all"), 
                                           x=trinity_aa[start:end], y=ensembl_aa[start:end]), ncol=5,byrow = TRUE))
        pbj$tick()
        
        
    }
    lcs = as.data.frame(lcs)
    return(lcs)
}

lcs_df_names = c("lcs_seq","a_start", "b_start","a_end","b_end")

if(do_lcs_frame %in%  c("p1", "p2", "p3")){
    aa1p = translate(DNAStringSet(trinity_translate$seq))
    aa2p = translate(DNAStringSet(str_sub(trinity_translate$seq,2,-1)) )
    aa3p = translate(DNAStringSet(str_sub(trinity_translate$seq,3,-1)) )
    trinity_translate_aas_pos = data.frame(trinity_id = trinity_translate$transcript_id, aa1p,aa2p,aa3p)
    
    
    mt = match(e2t$subject, trinity_translate_aas_pos$trinity_id)
    me = match(strv_split(e2t$query, "[.]", 1), strv_split(ensembl_peps$transcript, "[.]", 1))
    
    trinity_aa_m1p = trinity_translate_aas_pos$aa1p[mt]
    trinity_aa_m2p = trinity_translate_aas_pos$aa2p[mt]
    trinity_aa_m3p = trinity_translate_aas_pos$aa3p[mt]
    
    ensembl_aa = paste0(ensembl_peps$seq[me], "*")
}
if(do_lcs_frame ==  "p1"){
    lcs1p = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m1p, ensembl_aa = ensembl_aa);message("found LCS for frame 1")
    save(lcs1p, aa1p, trinity_translate_aas_pos, mt, me, ensembl_aa, trinity_aa_m1p, file=paste0("lcs_p1.rdata"))
}
if(do_lcs_frame ==  "p2"){
    lcs2p = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m2p, ensembl_aa = ensembl_aa);message("found LCS for frame 2")
    save(lcs2p, aa2p, trinity_translate_aas_pos, mt, me, ensembl_aa, trinity_aa_m2p, file=paste0("lcs_p2.rdata"))
}
if(do_lcs_frame ==  "p3"){
    lcs3p = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m3p, ensembl_aa = ensembl_aa);message("found LCS for frame 3")
    save(lcs3p, aa3p, trinity_translate_aas_pos, mt, me, ensembl_aa, trinity_aa_m3p, file=paste0("lcs_p3.rdata"))
}


if(do_lcs_frame %in%  c("n1", "n2", "n3")){
    aa1n = translate(reverseComplement(DNAStringSet(trinity_translate$seq)))
    aa2n = translate(reverseComplement(DNAStringSet(str_sub(trinity_translate$seq,1,-2))))
    aa3n = translate(reverseComplement(DNAStringSet(str_sub(trinity_translate$seq,1,-3))))
    trinity_translate_aas_neg = data.frame(trinity_id = trinity_translate$transcript_id, aa1n,aa2n,aa3n)
    
    mt = match(e2t$subject, trinity_translate_aas_neg$trinity_id)
    me = match(strv_split(e2t$query, "[.]", 1), strv_split(ensembl_peps$transcript, "[.]", 1))
    
    trinity_aa_m1n = trinity_translate_aas_neg$aa1n[mt]
    trinity_aa_m2n = trinity_translate_aas_neg$aa2n[mt]
    trinity_aa_m3n = trinity_translate_aas_neg$aa3n[mt]
    
    ensembl_aa = paste0(ensembl_peps$seq[me], "*")
}
if(do_lcs_frame ==  "n1"){
    lcs1n = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m1n, ensembl_aa = ensembl_aa);message("found LCS for frame 1")
    save(lcs1n, aa1n, trinity_translate_aas_neg, mt, me, ensembl_aa, trinity_aa_m1n, file=paste0("lcs_n1.rdata"))
}
if(do_lcs_frame ==  "n2"){
    lcs2n = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m2n, ensembl_aa = ensembl_aa);message("found LCS for frame 2")
    save(lcs2n, aa2n, trinity_translate_aas_neg, mt, me, ensembl_aa, trinity_aa_m2n, file=paste0("lcs_n2.rdata"))
}
if(do_lcs_frame ==  "n3"){
    lcs3n = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m3n, ensembl_aa = ensembl_aa);message("found LCS for frame 3")
    save(lcs3n, aa3n, trinity_translate_aas_neg, mt, me, ensembl_aa, trinity_aa_m3n, file=paste0("lcs_n3.rdata"))
}



