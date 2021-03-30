source("scripts/helper_functions.R")

######### go from cdna // cds differences

ensembl_upstream_stops = NULL
for(org in c('arabidopsis', 'zebrafish','human')){

    fasta_file = list.files(paste0("data/reference_annotations/", org, "/"), full.names = TRUE, pattern=".fa")
    cdna_file = fasta_file[grep('[.]cdna', fasta_file)]
    cds_file = fasta_file[grep('[.]cds', fasta_file)]
    
    cdna = halpme::read_fasta2df(cdna_file)
    cds = halpme::read_fasta2df(cds_file)
    
    cdna$transcript_id = strv_split(cdna$seq_id, "[ ]", 1)
    cds$transcript_id = strv_split(cds$seq_id, "[ ]", 1)
    
    m = match(cds$transcript_id, cdna$transcript_id)
    rm = which(cds$seq == cdna$seq[m])
    cds = cds[-rm,]
    m = match(cds$transcript_id, cdna$transcript_id)
    
    upstream_seq = as.character(mapply(function(x,y) str_sub(x, 1, (str_locate(x, y)[1] -1)), y=cds$seq, x=cdna$seq[m]))
    codon_change = nchar(upstream_seq) %% 3
    upstream_seq_adj = as.character(mapply(function(x,y) str_sub(y, x+1, -1), y=upstream_seq, x=codon_change))
    
    five_prime_utrs = data.frame(id = cds$transcript_id, seq = upstream_seq_adj)
    five_prime_utrs = five_prime_utrs[five_prime_utrs$seq != "",]
    five_prime_utrs = five_prime_utrs[which(!is.na(five_prime_utrs$seq)),]
    five_prime_utrs = five_prime_utrs[!grepl("N",five_prime_utrs$seq),]
    
    five_prime_utrs$aa_seq = as.character(translate(DNAStringSet(five_prime_utrs$seq)))
    
    #has_upstream_stop
    five_prime_utrs$aa_rev = stringi::stri_reverse(five_prime_utrs$aa_seq)
    stop_locs = str_locate(gsub("[*]","Z",five_prime_utrs$aa_rev),"Z")
    
    five_prime_utrs$upstream_stop = stop_locs[,1]
    five_prime_utrs$has_upstream_stop = ifelse(is.na(five_prime_utrs$upstream_stop), 0, 1)
    five_prime_utrs$upstream_stop_nt = five_prime_utrs$upstream_stop *3
    five_prime_utrs$upstream_stop_nt[five_prime_utrs$has_upstream_stop == 0] = nchar(five_prime_utrs$seq[five_prime_utrs$has_upstream_stop == 0])
    five_prime_utrs$upstream_stop[five_prime_utrs$has_upstream_stop == 0] = floor(five_prime_utrs$upstream_stop_nt[five_prime_utrs$has_upstream_stop == 0] /3)
    
    five_prime_utrs$utr5_len_nt = nchar(five_prime_utrs$seq)
    
    
    five_prime_utrs$seq = NULL
    five_prime_utrs$aa_rev = NULL
    
    write.table(five_prime_utrs, paste0("data/reference_annotations/", org, "/upstream_stops.txt"), quote=F, row.names = F, sep='\t')
    five_prime_utrs$org = org
    ensembl_upstream_stops = rbind(ensembl_upstream_stops, five_prime_utrs)
    
}

write.table(ensembl_upstream_stops, "data/ensembl_upstream_stops.txt", quote=F, row.names = F, sep='\t')





