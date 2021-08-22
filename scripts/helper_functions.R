library(progress)

combine_blastx_lines = function(blastx){

  ## combine multi blastx hits to single line
  qsub = paste0(blastx$query, blastx$subject)
  blastx_multi = blastx[qsub %in% qsub[which(duplicated(qsub))],]
  blastx_multi_summary = aggregate(p_ident ~ query + subject, blastx_multi, max)
  blastx_multi_summary = cbind(blastx_multi_summary, aggregate(length ~ query + subject, blastx_multi, sum)[,3])
  blastx_multi_summary = cbind(blastx_multi_summary, aggregate(mismatch ~ query + subject, blastx_multi, sum)[,3])
  blastx_multi_summary = cbind(blastx_multi_summary, aggregate(gapopen ~ query + subject, blastx_multi, sum)[,3])
  blastx_multi_summary$qstart = NA
  blastx_multi_summary$qend = NA
  blastx_multi_summary = cbind(blastx_multi_summary, aggregate(sstart ~ query + subject, blastx_multi, min)[,3])
  blastx_multi_summary = cbind(blastx_multi_summary, aggregate(send ~ query + subject, blastx_multi, max)[,3])
  blastx_multi_summary = cbind(blastx_multi_summary, aggregate(evalue ~ query + subject, blastx_multi, min)[,3])
  blastx_multi_summary = cbind(blastx_multi_summary, aggregate(bitscore ~ query + subject, blastx_multi, max)[,3])
  blastx_multi_summary$sframe = NA
  blastx_multi_summary = cbind(blastx_multi_summary, aggregate(qframe ~ query + subject, blastx_multi, function(x) paste0(sort(unique(x)), collapse = ";"))[,3])
  colnames(blastx_multi_summary) = colnames(blastx)

  blastx = rbind(blastx[!qsub %in% qsub[which(duplicated(qsub))],], blastx_multi_summary)

  return(blastx)
}

get_utr_coords_from_overlaps = function(utr_ol, utr_granges, ens_granges){

  # subtract CDS from 3- --wide UTR
  utr_matches = as.data.frame(utr_granges[utr_ol@to])
  ens_matches = as.data.frame(ens_granges[utr_ol@from])

  # get 'UTR' coords by adding/subtracting from the end of the matched gene body
  left = utr_matches
  left$end = ens_matches$start - 1
  left = left[-which(left$start == left$end + 1),]
  left$width = NULL
  left = distinct(left)
  # Take only 1 UTR for each gene start/end -- the longest
  left = aggregate(start ~ seqnames+end+strand, left, min)[,c(1,4,2,3)] %>% arrange(seqnames, start, end)
  left_gr = GRanges(seqnames = left$seqnames, strand = left$strand, ranges=IRanges(start = left$end+1, end = left$end+1))
  left_ol = findOverlaps(left_gr, ens_granges, type = "start")
  left = left[left_ol@from,]
  left$gene_id = ens_granges$transcript[left_ol@to]

  right = utr_matches
  right$start = ens_matches$end + 1
  right = right[-which(right$start == right$end + 1),]
  right$width = NULL
  right = distinct(right)
  # Take only 1 UTR for each gene start/end -- the longest
  right = aggregate(end ~ seqnames+start+strand, right, max)[,c(1,2,4,3)] %>% arrange(seqnames, start, end)
  right_gr = GRanges(seqnames = right$seqnames, strand = right$strand, ranges=IRanges(start = right$start-1, end = right$start-1))
  right_ol = findOverlaps(right_gr, ens_granges, type = "end")
  right = right[right_ol@from,]
  right$gene_id = ens_granges$transcript[right_ol@to]

  # convert left/right to 3'/5'
  utr3 = rbind(left[left$strand == "+",], right[right$strand == "-",])
  utr5 = rbind(left[left$strand == "-",], right[right$strand == "=",])
  utr3$type = "three_prime_UTR"
  utr5$type = "five_prime_UTR"

  utrs = rbind(utr5, utr3)
  utrs$transcript = paste0(utrs$gene_id, "_mRNA")

  return(utrs)
}

get_utr_cds_summary = function(ensembl_gff_utr3, ensembl_gff_utr5, ensembl_gff_cds, ens_granges, utr_granges){
  ##### convert ensembl gffs to CDS locations relative to tx start/end
  #ensembl_gff_utr3 = utrs[utrs$type == "three_prime_UTR",]
  #ensembl_gff_utr5 = utrs[utrs$type == "three_prime_UTR",]

  ensembl_gff_utr5$width = (ensembl_gff_utr5$end - ensembl_gff_utr5$start) +1
  ensembl_gff_utr3$width = (ensembl_gff_utr3$end - ensembl_gff_utr3$start) +1

  utr5_lens = aggregate(width ~ transcript, ensembl_gff_utr5, sum)
  utr3_lens = aggregate(width ~ transcript, ensembl_gff_utr3, sum)

  ens_cds_summary = aggregate(width ~ transcript, ensembl_gff_cds, sum)
  ens_cds_summary$utr5_len = NA
  ens_cds_summary$utr3_len =NA

  utr_ol_evidence = findOverlaps(ens_granges, utr_granges)
  genes_with_utr_evidence = unique(ens_granges$transcript[utr_ol_evidence@from])
  ens_cds_summary$utr5_len[ens_cds_summary$transcript %in% paste0(genes_with_utr_evidence,"_mRNA")] = 0
  ens_cds_summary$utr3_len[ens_cds_summary$transcript %in% paste0(genes_with_utr_evidence,"_mRNA")] = 0

  ens_cds_summary$utr5_len[which(!is.na(match(ens_cds_summary$transcript, utr5_lens$transcript)))] =
    utr5_lens$width[match(ens_cds_summary$transcript, utr5_lens$transcript)][which(!is.na(match(ens_cds_summary$transcript, utr5_lens$transcript)))]
  ens_cds_summary$utr3_len[which(!is.na(match(ens_cds_summary$transcript, utr3_lens$transcript)))] =
    utr3_lens$width[match(ens_cds_summary$transcript, utr3_lens$transcript)][which(!is.na(match(ens_cds_summary$transcript, utr3_lens$transcript)))]
  return(ens_cds_summary)
}


#a = "SDGFADSAGDGAXYZXYZXYZDGHAJFD"
#b = "XYZXYZXYZDGHJFKGDSHJFK"

longest_common_substring = function(a,b, return='seq'){
  if(is.na(a) | is.na(b) | nchar(a) < 4 | nchar(b) < 4){

    lcs_df = data.frame(lcs_seq = NA,
                        a_start = NA,
                        b_start = NA,
                        a_end = NA,
                        b_end = NA)
    if(return == 'seq'){
      return(NA)
    }else{
      return(lcs_df)
    }

  }else{
    A <- strsplit(a, "")[[1]]
    B <- strsplit(b, "")[[1]]

    L <- matrix(0, length(A), length(B))
    ones <- which(outer(A, B, "=="), arr.ind = TRUE)
    ones <- ones[order(ones[, 1]), ]

    onesd = as.data.frame(ones)
    onesd$diff = onesd$col - onesd$row
    onesd = arrange(onesd, desc(diff), row, col)
    onesd$same_as_prev = c(0, onesd$col[-1] - onesd$col[-nrow(onesd)]) ==1 & c(0, onesd$row[-1] - onesd$row[-nrow(onesd)]) == 1
    len = c(which(onesd$same_as_prev == F)[-1] -1, nrow(onesd)) - which(onesd$same_as_prev == F) +1

    onesd = onesd[which(onesd$same_as_prev == F),c(1,2)]
    onesd$len = len
    onesd = onesd[which.max(onesd$len),]

    lcs_df = data.frame(lcs_seq = paste0(A[onesd$row:(onesd$row+onesd$len-1)], collapse =""),
                        a_start = onesd$row[1],
                        b_start = onesd$col[1],
                        a_end = onesd$row[1] + onesd$len[1] -1,
                        b_end = onesd$col[1] + onesd$len[1] -1)

    if(return == 'seq'){
      return(lcs_df$lcs_seq)
    }else{
      return(lcs_df)
    }
  }
}


get_lcs_trinity_ensembl = function(trinity_translate, ensembl_peps, e2t, chunk_size = NA, strand = "+"){

  make_lcs_df_chunk = function(chunk_size, trinity_aa, ensembl_aa){

    chunks = ceiling(length(trinity_aa) / chunk_size)
    lcs = NULL
    pbj = progress_bar$new(total = chunks, format = "[:bar] :current/:total (:percent) time: :elapsedfull eta: :eta")
    for(i in 1:chunks){

      start = i*chunk_size - chunk_size + 1
      end = min(length(trinity_aa), start+chunk_size-1)
      lcs = rbind(lcs, matrix(mapply(function(x,y) longest_common_substring(x, y, return="all"),
                                       x=trinity_aa[start:end], y=ensembl_aa[start:end]), ncol=5,byrow = TRUE))
      pbj$tick()


    }
    lcs = as.data.frame(lcs)
    return(lcs)
  }

  lcs_df_names = c("lcs_seq","a_start", "b_start","a_end","b_end")
  if(strand == "+" | strand == "both"){
    aa1p = translate(DNAStringSet(trinity_translate$seq))
    aa2p = translate(DNAStringSet(str_sub(trinity_translate$seq,2,-1)) )
    aa3p = translate(DNAStringSet(str_sub(trinity_translate$seq,3,-1)) )
    trinity_translate_aas_pos = data.frame(trinity_id = trinity_translate$transcript_id, aa1p,aa2p,aa3p)


    mt = match(e2t$subject, trinity_translate_aas_pos$trinity_id)
    me = match(e2t$query, ensembl_peps$transcript)

    trinity_aa_m1p = trinity_translate_aas_pos$aa1p[mt]
    trinity_aa_m2p = trinity_translate_aas_pos$aa2p[mt]
    trinity_aa_m3p = trinity_translate_aas_pos$aa3p[mt]

    ensembl_aa = paste0(ensembl_peps$seq[me], "*")
    if(!is.na(chunk_size)){
      lcs1p = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m1p, ensembl_aa = ensembl_aa);message("found LCS for frame 1")
      lcs2p = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m2p, ensembl_aa = ensembl_aa);message("found LCS for frame 2")
      lcs3p = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m3p, ensembl_aa = ensembl_aa);message("found LCS for frame 3")
    }else{
      lcs1p = matrix(mapply(function(x,y) longest_common_substring(x, y, return="all"), x=trinity_aa_m1p, y=ensembl_aa), ncol=5,byrow = TRUE)
      message("found LCS for frame 1")
      lcs2p = matrix(mapply(function(x,y) longest_common_substring(x, y, return="all"), x=trinity_aa_m2p, y=ensembl_aa), ncol=5,byrow = TRUE)
      message("found LCS for frame 2")
      lcs3p = matrix(mapply(function(x,y) longest_common_substring(x, y, return="all"), x=trinity_aa_m3p, y=ensembl_aa), ncol=5,byrow = TRUE)
      message("found LCS for frame 3")
    }

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
  }

  if(strand == "-" | strand == "both"){
    aa1n = translate(reverseComplement(DNAStringSet(trinity_translate$seq)))
    aa2n = translate(reverseComplement(DNAStringSet(str_sub(trinity_translate$seq,2,-1))))
    aa3n = translate(reverseComplement(DNAStringSet(str_sub(trinity_translate$seq,3,-1))))
    trinity_translate_aas_neg = data.frame(trinity_id = trinity_translate$transcript_id, aa1n,aa2n,aa3n)

    mt = match(e2t$subject, trinity_translate_aas_neg$trinity_id)
    me = match(e2t$query, ensembl_peps$transcript)

    trinity_aa_m1n = trinity_translate_aas_neg$aa1n[mt]
    trinity_aa_m2n = trinity_translate_aas_neg$aa2n[mt]
    trinity_aa_m3n = trinity_translate_aas_neg$aa3n[mt]

    ensembl_aa = paste0(ensembl_peps$seq[me], "*")
    if(!is.na(chunk_size)){
      lcs1n = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m1n, ensembl_aa = ensembl_aa);message("found LCS for frame 1")
      lcs2n = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m2n, ensembl_aa = ensembl_aa);message("found LCS for frame 2")
      lcs3n = make_lcs_df_chunk(chunk_size=chunk_size, trinity_aa = trinity_aa_m3n, ensembl_aa = ensembl_aa);message("found LCS for frame 3")
    }else{
      lcs1n = matrix(mapply(function(x,y) longest_common_substring(x, y, return="all"), x=trinity_aa_m1n, y=ensembl_aa), ncol=5,byrow = TRUE)
      message("found LCS for frame 1")
      lcs2n = matrix(mapply(function(x,y) longest_common_substring(x, y, return="all"), x=trinity_aa_m2n, y=ensembl_aa), ncol=5,byrow = TRUE)
      message("found LCS for frame 2")
      lcs3n = matrix(mapply(function(x,y) longest_common_substring(x, y, return="all"), x=trinity_aa_m3n, y=ensembl_aa), ncol=5,byrow = TRUE)
      message("found LCS for frame 3")
    }

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
  }




  if(strand == "+"){
    lcs_df = lcs_dfp
  }else if(strand == "-"){
    lcs_df = lcs_dfn
  }else{

    #lcs_dfp$lcs_len = nchar(lcs_dfp$lcs_seq)
    #lcs_dfn$lcs_len = nchar(lcs_dfn$lcs_seq)

    pn = mapply(function(x,y) which.max(c(x,y))[1], x=lcs_dfp$lcs_len, y=lcs_dfn$lcs_len)
    lcs_df = lcs_dfp
    lcs_df[which(pn == 2),] = lcs_dfn[which(pn == 2),]

  }


  lcs_df$lcs_seq = as.character(unlist(lcs_df$lcs_seq))
  for(i in 7:10){
    lcs_df[,i] = as.numeric(unlist(lcs_df[,i]))
  }


  lcs_df$first_AA = ifelse(str_sub(lcs_df$lcs_seq,1,1)=="M", "MET", "ALT")
  lcs_df$last_AA = ifelse(str_sub(lcs_df$lcs_seq,-1,-1)=="*", "STOP", "ALT")
  lcs_df$lcs_len = ifelse(lcs_df$last_AA == "STOP", lcs_df$lcs_len-1, lcs_df$lcs_len)

  return(lcs_df[,-c(1:3)])

}

read_lst= function(filename){
    lst = fread(cmd=paste0("cat ", filename, ' | grep -v "Predicted genes" | grep -v "Model information" | grep -v "LeftEnd" | grep -v "#" | grep . '), skip = 6, header=F, sep='\n')
    fa_def_lines = grep("FASTA definition line", lst$V1)
    fa_def_len = fa_def_lines[-1] - fa_def_lines[-length(fa_def_lines)]

    rm = fa_def_lines[fa_def_len == 1]

    lst_wide = NULL
    for(i in 2:max(fa_def_len)){
        if(i == 2){
            lst_wide = rbind(lst_wide, data.frame(
                seq_id = lst[fa_def_lines[fa_def_len == i],],
                seq_feats = lst[fa_def_lines[fa_def_len == i]+1,]
            ))
        }else{
            seq_id = lst[rep(fa_def_lines[which(fa_def_len == i)], each=(i-1)),]
            seq_feats = lst[rep(fa_def_lines[which(fa_def_len == i)], each=(i-1))+
                                rep(c(1:(i-1)),length(which(fa_def_len == i))),]
            lst_wide = rbind(lst_wide, data.frame(seq_id, seq_feats))
        }
    }

    #  GeneNumber    Strand    LeftEnd    RightEnd   GeneLength    Class
    lst_wide$V1.1 = str_squish(lst_wide$V1.1)
    n_spaces = unlist(lapply(str_split(lst_wide$V1.1, "[ ]"), function(x) length(x)))
    lst_wide = lst_wide[n_spaces==6,]

    lst_wide$gene_number = strv_split(lst_wide$V1.1, "[ ]", 1)
    lst_wide$strand = strv_split(lst_wide$V1.1, "[ ]", 2)
    lst_wide$leftend = strv_split(lst_wide$V1.1, "[ ]", 3)

    lst_wide$rightend = strv_split(lst_wide$V1.1, "[ ]", 4)

    lst_wide$gene_length = strv_split(lst_wide$V1.1, "[ ]", 5)
    lst_wide$class = strv_split(lst_wide$V1.1, "[ ]", 6)

    lst_wide$V1 = gsub("FASTA definition line: ", "", lst_wide$V1)
    colnames(lst_wide)[1] = "seq_id"
    lst_wide$V1.1=NULL
    lst_wide$orf_class = "complete"
    lst_wide$orf_class[grep("<", lst_wide$leftend)] = "incomplete_5prime"
    lst_wide$orf_class[grep(">", lst_wide$rightend)] = "incomplete_3prime"
    lst_wide$orf_class[grepl(">", lst_wide$rightend) & grepl("<", lst_wide$leftend)] = "incomplete"
    return(lst_wide)
}
