#check splice site alignment

library(data.table)
library(tidyverse)
library(stringr)
library(halpme) # github: betsig/halpme
library(GenomicRanges)
library(stringdist)
library(Biostrings)

options(stringsAsFactors = F)

org = 'human'
exonsToIntrons <- function(exons) {
    exons_df <- as.data.frame(exons)
    exons_df <- exons_df[, c("seqid", "start", "end", "strand", "transcript_id", "exon_number")]
    exons_df <- exons_df[order(exons_df$transcript_id, exons_df$start, exons_df$end), ]
    exons_df$intron_start <- exons_df$end
    exons_df$intron_end <- dplyr::lead(exons_df$start)
    
    rm <- which(dplyr::lead(exons_df$transcript_id) != exons_df$transcript_id)
    if(length(rm) > 0){
        exons_df <- exons_df[-rm, ]
    }
    exons_df <- exons_df[-nrow(exons_df), ]
    min_exon_n <- aggregate(exon_number ~ transcript_id, exons_df, min)
    if (!all(min_exon_n$exon_number == 1)) {
        exons_df$exon_number[exons_df$strand == "-"] <- as.numeric(exons_df$exon_number[exons_df$strand == "-"]) - 1
    }
    
    introns <- GRanges(
        seqnames = exons_df$seqid, ranges = IRanges(start = exons_df$intron_start, end = exons_df$intron_end),
        strand = exons_df$strand, transcript_id = exons_df$transcript_id, exon_number = exons_df$exon_number
    )
    m <- match(introns$transcript_id, exons$transcript_id)
    introns$gene_id <- exons$gene_id[m]
    introns$gene_name <- exons$gene_name[m]
    
    #mcols(introns) <- DataFrame(dplyr::left_join(as.data.frame(mcols(introns)), exons))
    return(introns)
}





tt_exp_all = NULL
for(org in c("human", "arabidopsis", "zebrafish")){

    load(paste0("data/", org, "/processed_ensembl.Rdata"))
    
    tt = fread(paste0('data/', org, "/trinity_translated_both_strands.txt"), data.table = F)
    trinity_seqs = fread(paste0("data/", org, "/trinity_with_kallisto.txt"), data.table = F)
    colnames(trinity_seqs) = c("full_id", "subject", "trinity_gene", "trinity_length_nt", "read_counts", "reads_per_kb")
    trinity_seqs$base_class = "no_hits"
    tt = plyr::rbind.fill(tt, trinity_seqs[which(!(trinity_seqs$subject %in% tt$subject)),-1])
    tt$org = org
    #trinity_translated = plyr::rbind.fill(trinity_translated, tt)
    fasta_file=ifelse(file.exists(paste0('data/', org, "/Trinity.fasta.gz")),paste0('data/', org, "/Trinity.fasta.gz"), paste0('data/', org, "/Trinity.fasta"))
    trinity_fasta = halpme::read_fasta2df(fasta_file)
    trinity_fasta$trinity_id = strv_split(trinity_fasta$seq_id, "[ ]", 1)
    
    tt_expressed = tt[which(tt$read_counts >=100),]
    #tt_expressed$ens_gene = ensembl_gff$gene_id[match(strv_split(tt_expressed$query_gene, "[.]", 1), ensembl_gff$transcript_id)]
    
    ## try cdna headers... 
    
    ens_cdna_headers = halpme::read_fasta2df(list.files(paste0('data/reference_annotations/',org), pattern="cdna", full.names = T)[1])[,1]
    ens_cdna = data.frame(transcript_id = strv_split(ens_cdna_headers, "[ ]" , 1),
                          #gene_id = strv_split2(ens_cdna_headers, "gene[:]", "[ ]"), 
                          loc = strv_split(ens_cdna_headers, "[ ]", 3))
    if(all(grepl("gene[:]", ens_cdna_headers[1:10]))){
        ens_cdna$gene_id = strv_split2(ens_cdna_headers, "gene[:]", "[ ]")
    }else{
        ens_cdna$gene_id = strv_split(ens_cdna_headers, "[|]", 2)
        ens_cdna$loc = strv_split(ens_cdna_headers, "[ ][|][ ]", 3)
    }
    if(org == "zebrafish" | org == "human"){
        ens_cdna$chrom=strv_split(ens_cdna$loc, "[:]", 3)
    }else{
        ens_cdna$chrom = strv_split(ens_cdna$loc, ":", 1)
    }
    if(org == "arabidopsis"){
        ens_cdna$gene_id = ensembl_gff$gene_id[match(ens_cdna$transcript_id, ensembl_gff$transcript_id)]
    }
    
    m1 = match(tt_expressed$query_gene, ens_cdna$transcript_id)
    tt_expressed$ens_gene = ens_cdna$gene_id[m1]
    tt_expressed$ens_transcript = ens_cdna$transcript_id[m1]
    tt_expressed$ens_chr = ens_cdna$chrom[m1]
    tt_expressed = tt_expressed[which(tt_expressed$ens_chr %in% c(1:22, "MT", "X", "Y", "C", "M", paste0("Chr", c(1:22, "MT", "X", "Y", "C", "M")))),]

    gmap = fread(paste0("data/", org, "/Trinity.gmap.gff"), data.table = F)
    colnames(gmap) = c("seqname", "source", "type", "start", "end", "score", "strand", "frame", "attributes")
    gmap$trinity_id = strv_split2(gmap$attributes, "Name[=]", ";")
    gmap$trinity_map_id = strv_split2(gmap$attributes, "ID[=]", ";")
    gmap$trinity_start = strv_split(strv_split2(gmap$attributes, "Target[=]", ";"), "[ ]", 2)
    gmap$trinity_end = strv_split(strv_split2(gmap$attributes, "Target[=]", ";"), "[ ]", 3)
    gmap$path_number = strv_split(gmap$trinity_map_id, "path", 2)
    
    
    if(org == "arabidopsis"){
        gmap$seqname[gmap$seqname == "chloroplast"] = "C"
        gmap$seqname[gmap$seqname == "mitochondria"] = "M"
        gmap$seqname = paste0("Chr", gmap$seqname)
    }

    gene_coords = aggregate(start ~ gmap$trinity_map_id, gmap, min)
    gene_coords$end = aggregate(end ~ gmap$trinity_map_id, gmap, max)[,2]
    colnames(gene_coords)[1] = "trinity_map_id"
    gene_coords = cbind(gene_coords, gmap[match(gene_coords$trinity_map_id, gmap$trinity_map_id), -c(4,5,11:14)])
    gene_coords$blast_match_transcript = strv_split(tt_expressed$query_gene[match(gene_coords$trinity_id, tt_expressed$subject)], "[.]", 1)
    gene_coords$blast_match_gene = strv_split(tt_expressed$ens_gene[match(gene_coords$trinity_id, tt_expressed$subject)], "[.]", 1)
    
    
    gmap_granges = GRanges(seqnames = gene_coords$seqname, ranges=IRanges(start=gene_coords$start, end=gene_coords$end), strand = gene_coords$strand, 
                           trinity_map_id = gene_coords$trinity_map_id, trinity_id = gene_coords$trinity_id, blast_gene = gene_coords$blast_match_gene)
    
    ensembl_granges = GRanges(seqnames = ensembl_gff_gene$seqid, ranges = IRanges(ensembl_gff_gene$start, ensembl_gff_gene$end), strand = ensembl_gff_gene$strand, 
                              gene_id = ensembl_gff_gene$gene_id)
    
    
    # check overlaps with gmapped coords and ref genes
    gmap2ens = as.data.frame(findOverlaps(gmap_granges, ensembl_granges, ignore.strand = T))

    gmap2ens$gmap_blast_hit = gmap_granges$blast_gene[gmap2ens$queryHits]
    gmap2ens$ens_gene = ensembl_granges$gene_id[gmap2ens$subjectHits]
    gmap2ens$trinity_path_id = gene_coords$trinity_map_id[gmap2ens$queryHits]
    gmap2ens$trinity_id = gene_coords$trinity_id[gmap2ens$queryHits]
    
    gmap2ens.correct = gmap2ens[which(gmap2ens$gmap_blast_hit == gmap2ens$ens_gene),]
    
    table(tt$subject %in% gene_coords$trinity_id)
    
    # for each gmapped transcript, check if it maps to the right gene (the best blast hit)
    gene_coords$GMAPPED_CORRECT = ifelse(gene_coords$trinity_map_id %in% gmap2ens.correct$trinity_path_id, T, F)
    # some gmapped transcripts may multi map (especialy in the case of antisense txs)
    gene_coords = arrange(gene_coords, desc(GMAPPED_CORRECT))
    tt_expressed$GMAPPED_CORRECT = ifelse(tt_expressed$subject %in% gene_coords$trinity_id[gene_coords$GMAPPED_CORRECT==T], T,F)
    
    table(tt_expressed$GMAPPED_CORRECT, tt_expressed$base_class)
    
    #### find splice sites and check overlaps
    
    # transcripts where tx maps to correct gene
    test = tt_expressed[tt_expressed$GMAPPED_CORRECT == T,]
    # gmap coords for each correctly mapped tx    
    gmap_match = gmap[gmap$trinity_map_id %in% gmap2ens.correct[match(test$subject, gmap2ens.correct$trinity_id),]$trinity_path_id,]
    
    # add in exon numbers
    gmap_match$exon_number = NA
    wna = which(is.na(gmap_match$exon_number))
    exon_num=1
    while(length(wna)>0){
        gmap_match$exon_number[wna][which(!(duplicated(gmap_match$trinity_map_id[wna])))] = exon_num
        wna = which(is.na(gmap_match$exon_number))
        exon_num=exon_num+1
    }
    
    colnames(gmap_match)[1] = "seqid"    
    colnames(gmap_match)[c(10,11)] = c("gene_id", "transcript_id")
        
    #exons = ensembl_gff[ensembl_gff$transcript_id %in% strv_split(test$query, "[.]", 1),]
    exons = ensembl_gff[ensembl_gff$gene_id %in% strv_split(test$ens_gene, "[.]", 1),]
    if(all(grepl("exon_number", exons$attributes))){
        exons$exon_number = strv_split2(exons$attributes, 'exon_number "', '";')
    }else{
        exons$exon_number = strv_split2(exons$attributes, 'exon[:]', ";")
        exons = separate_rows(exons, transcript_id, sep=',')
        exons$gene_name = exons$gene_id
    }
    # convert exons 2 splice junctions
    introns = exonsToIntrons(exons)
    gmap_introns = exonsToIntrons(gmap_match)
            
    # find identical splice site pairs
    ol.introns = findOverlaps(gmap_introns, introns, type="equal", ignore.strand=T)
    ol.introns = as.data.frame(ol.introns)
    
    # double check that overlapping gene ids match
    ol.introns$from_trinity = gmap_introns$gene_id[ol.introns$queryHits]
        #ol.introns$to_ensembl = introns$transcript_id[ol.introns$subjectHits]
    ol.introns$to_ensembl = introns$gene_id[ol.introns$subjectHits]
        #ol.introns$correct_ens = strv_split(test$query[match(ol.introns$from_trinity, test$subject)], "[.]", 1)
    ol.introns$correct_ens = strv_split(test$ens_gene[match(ol.introns$from_trinity, test$subject)], "[.]", 1)
    ol.introns.same = ol.introns[ol.introns$to_ensembl == ol.introns$correct_ens,]
        
    n_intron_correct = as.data.frame(table(ol.introns.same$from_trinity[!duplicated(paste0(ol.introns.same$queryHits, "_", ol.introns.same$from_trinity))]))
    total_introns = as.data.frame(table(gmap_introns$transcript_id))
    total_introns$trinity_id = strv_split(total_introns$Var1, "[.]", 1)
    total_introns$n_correct = n_intron_correct$Freq[match(total_introns$trinity_id,n_intron_correct$Var1)]
    total_introns$n_correct[which(is.na(total_introns$n_correct))] = 0
    total_introns$percent_correct = total_introns$n_correct / total_introns$Freq
    n_intron_correct$percent_correct = n_intron_correct$Freq/total_introns$Freq[match(n_intron_correct$Var1, total_introns$trinity_id)]
        
    tt_expressed$n_intron_correct = total_introns$n_correct[match(tt_expressed$subject, strv_split(total_introns$Var1, "[.]" ,1))]
    tt_expressed$percent_intron_correct = total_introns$percent_correct[match(tt_expressed$subject, strv_split(total_introns$Var1, "[.]" ,1))]
    
    tt_exp_all = rbind(tt_exp_all, tt_expressed)
}
tt_exp_all$lcs_seq=NULL
tt_exp_all$trinity_orf=NULL
write.table(tt_exp_all, "data/splice_site_checks2.txt", row.names = F, sep='\t')

ggplot(tt_exp_all, aes(x=org, y=percent_intron_correct, fill=match_strand)) + 
    geom_violin() + ggeasy::easy_rotate_x_labels()

ggplot(tt_expressed, aes(fill=base_class, x=percent_intron_correct)) + geom_density(alpha=0.2) + facet_wrap(~match_strand) + ggeasy::easy_rotate_x_labels()

aggregate(percent_intron_correct~base_class+match_strand, tt_expressed, count)

table(tt_expressed$match_strand, tt_expressed$percent_intron_correct ==0) / rowSums(table(tt_expressed$match_strand, tt_expressed$percent_intron_correct ==0))

x = table(tt_exp_all$match_strand, tt_exp_all$percent_intron_correct ==0, tt_exp_all$org) %>% as.data.frame() %>% reshape(idvar = c("Var1", "Var3"), timevar = "Var2", direction = "wide")
x$percent_with_junc = x$Freq.FALSE/rowSums(x[,c(3,4)])


total_introns$match_strand = tt_expressed$match_strand[match(total_introns$trinity_id, tt_expressed$subject)]



### double check a false negative strand tx works
fake = 'TRINITY_DN2759_c2_g1_i4'
tt.insert = tt_expressed[tt_expressed$subject == fake,]
tt.insert$match_strand == "-"
tt.insert$subject = "TEST"

gmap.insert = gmap[gmap$trinity_id == fake,]
gmap.insert$strand == "-"
gmap.insert$trinity_id = "TEST"



tt_expressed %>% filter(GMAPPED_CORRECT == F & base_class == "complete") %>% head(1) -> incor


gmap2ens[gmap2ens$trinity_id==incor$subject,]

ol = findOverlaps(gmap_granges[grep(incor$subject, gmap_granges$trinity_map_id)], ensembl_granges, ignore.strand=T)
gmap_granges[grep(incor$subject, gmap_granges$trinity_map_id)][ol@from]
ensembl_granges[which(ensembl_granges$gene_id %in% gmap_granges[grep(incor$subject, gmap_granges$trinity_map_id)][ol@from]$blast_gene)]


