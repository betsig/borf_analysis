# strandedness plots

library(data.table)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(stringdist)
library(Biostrings)
if (!require('halpme')) devtools::install_github("betsig/halpme")
library(halpme)
if (!require('ggeasy')) devtools::install_github("jonocarroll/ggeasy")
library(ggeasy)

# main figure theme
theme_figure <- theme_bw()+ theme(text=element_text(size=8),legend.key.size=unit(0.1, "inches"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_rect(colour = "black", size=0.5, fill=NA),
                                  axis.line.x = element_blank(),
                                  axis.line.y = element_blank(),
                                  axis.ticks = element_line(size=0.25),
                                  axis.ticks.length = unit(0.05, "cm"),
                                  panel.background = element_blank(),
                                  plot.title = element_text(hjust = 0.5),
                                  strip.background = element_blank(),
                                  strip.text = element_text(size=8, face = "bold"))

# function to extract the legend from a ggplot object
ggplot_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

options(stringsAsFactors = F)

at_color_pals = read.delim("at_colour_pal.txt")
base_class_pal_muted = c("#fb8072","#fccde5","#ffffb3","#d9d9d9","#ccebc5","#bc80bd","#ffed6f","#80b1d3","#b3de69")
base_class_pal = c("#e31a1c","#cab2d6","#fdbf6f","#a6cee3","#b2df8a","#6a3d9a","#ff7f00","#1f78b4","#33a02c")
species_pal = c("#36B37E","#00B8D9","#FF5630","#FFAB00")

inten = 300
base_class_pal2 = c(
  at_color_pals$hex_value[at_color_pals$colour == "N" & at_color_pals$intensity ==inten+100],
  at_color_pals$hex_value[at_color_pals$colour == "P" & at_color_pals$intensity ==inten],
  at_color_pals$hex_value[at_color_pals$colour == "G" & at_color_pals$intensity ==inten-100],
  at_color_pals$hex_value[at_color_pals$colour == "B" & at_color_pals$intensity ==inten],
  at_color_pals$hex_value[at_color_pals$colour == "T" & at_color_pals$intensity ==inten],
  at_color_pals$hex_value[at_color_pals$colour == "G" & at_color_pals$intensity ==inten],
  at_color_pals$hex_value[at_color_pals$colour == "Y" & at_color_pals$intensity ==inten],
  at_color_pals$hex_value[at_color_pals$colour == "R" & at_color_pals$intensity ==inten])




######### Read in annotations from standard and alternative datasets ##########


alt_orgs = c('human', 'human', 'yeast', 'yeast', "human", "yeast","zebrafish", "arabidopsis")
alt_versions = c("cds", "pc", "2sample", "3sample", "","","","")


# need to read in both the strand-agnostic (bs/both strands) and the positve strand only annotations

trinity_translated_bs = NULL
for(i in 1:length(alt_orgs)){
    
    org = alt_orgs[i]
    version = alt_versions[i]
    org_version = ifelse(version=="", org, paste(org, version, sep="_"))
    
    tt = fread(paste0("data/", org_version, "/trinity_translated_both_strands.txt"), data.table = F)
    tt$org = org_version
    if(org == "human"){
        # human uses the same initial Trinity assembly for the cds/pc only blast-matching annotations
        coding_potential = fread(paste0("data/", org, "/rnasamba/rnasamba_partial.tsv"), data.table = F)
    }else{
        coding_potential = fread(paste0("data/", org_version, "/rnasamba/rnasamba_partial.tsv"), data.table = F)
    }
    tt$coding_potential = coding_potential$coding_score[match(tt$subject, halpme::strv_split(coding_potential$sequence_name, "[ ]", 1))]
    trinity_translated_bs = plyr::rbind.fill(trinity_translated_bs, tt)
}

trinity_translated = NULL
for(i in 1:length(alt_orgs)){
    
    org = alt_orgs[i]
    version = alt_versions[i]
    org_version = ifelse(version=="", org, paste(org, version, sep="_"))
    
    
    tt = fread(paste0("data/", org_version, "/trinity_translated.txt"), data.table = F)
    tt$org = org_version
    if(org == "human"){
        # human uses the same initial Trinity assembly for the cds/pc only blast-matching annotations
        coding_potential = fread(paste0("data/", org, "/rnasamba/rnasamba_partial.tsv"), data.table = F)
    }else{
        coding_potential = fread(paste0("data/", org_version, "/rnasamba/rnasamba_partial.tsv"), data.table = F)
    }
    tt$coding_potential = coding_potential$coding_score[match(tt$subject, halpme::strv_split(coding_potential$sequence_name, "[ ]", 1))]
    trinity_translated = plyr::rbind.fill(trinity_translated, tt)
}

trinity_translated_bs$base_class_factor = factor(trinity_translated_bs$base_class, levels = rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage")))
trinity_translated$base_class_factor = factor(trinity_translated$base_class, levels = rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage")))

add_organism_col = function(df){
    
    df$organism = "h. sapiens"
    df$organism[df$org == "zebrafish"] = "d. rerio"
    df$organism[df$org == "arabidopsis"] = "a. thaliana"
    df$organism[df$org == "yeast"] = "s. cerevisiae"
    df$organism[df$org == "yeast_2sample"] = "s. cerevisiae - 2 samples (116M reads)"
    df$organism[df$org == "yeast_3sample"] = "s. cerevisiae - 3 samples (171M reads)"
    df$organism[df$org == "yeast_dataset2"] = "s. cerevisiae - alt. dataset 2"
    df$organism[df$org == "yeast_dataset3"] = "s. cerevisiae - alt. dataset 3"
    
    df$organism[df$org == "human_cds"] = "h. sapiens - CDS only"
    df$organism[df$org == "human_pc"] = "h. sapiens - protein_coding only"
    
    return(df)
    
}

trinity_translated_bs = add_organism_col(trinity_translated_bs)
trinity_translated = add_organism_col(trinity_translated)

trinity_translated_bs$base_class_f2 = factor(gsub("_", " ", trinity_translated_bs$base_class), levels = c(levels = gsub("_"," ",rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage")))))
trinity_translated_bs$match_strand = factor(trinity_translated_bs$match_strand, levels = c("+",  "-"))

trinity_translated$base_class_f2 = factor(gsub("_", " ", trinity_translated$base_class), levels = c(levels = gsub("_"," ",rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage")))))
trinity_translated$match_strand = factor(trinity_translated$match_strand, levels = c("+",  "-"))


orgs = c("yeast", "human", "zebrafish", "arabidopsis")
org_scis=c("s. cerevisiae", "h. sapiens", "d. rerio", "a. thaliana")

all_tx_lengths = NULL
negative_annots_with_dists = NULL
all_gene_dists = NULL
gene_numbers = NULL

# for each species, get:
# number of reference genes, distance between reference genes
# distances from -ve blast match transcripts to antisense reference CDS

for(org in c("yeast", "human", "zebrafish", "arabidopsis")){
  org_sci = org_scis[which(orgs==org)]
  load(paste0("data/",org,"/processed_ensembl.Rdata"))
  if(org == "arabidopsis"){
    ## IF ARABIDOPSIS
    ensembl_gff = separate_rows(ensembl_gff, transcript_id)
  }
  
  # make 'genes' with CDS only 
  gene_cds = aggregate(start ~ gene_id, ensembl_gff_cds, min)
  gene_cds$end = aggregate(end ~ gene_id, ensembl_gff_cds, max)[,2]
  gene_cds = cbind(gene_cds, ensembl_gff_cds[match(gene_cds$gene_id, ensembl_gff_cds$gene_id),c('seqid','source', 'type','score','strand', 'phase')])
  genes = GenomicRanges::makeGRangesFromDataFrame(gene_cds)
  genes$gene_id = gene_cds$gene_id
  
  ## DISTANCES BETWEEN CDS's IN ANTISENSE GENES
  genes.pos = genes[strand(genes) == "+"]
  genes.neg = genes[strand(genes) == "-"]
  between_gene_dists.1 = GenomicRanges::distanceToNearest(genes.neg, genes.pos, ignore.strand = T)
  between_gene_dists.2 = GenomicRanges::distanceToNearest(genes.pos, genes.neg, ignore.strand = T)
  gene_dists = c(as.data.frame(between_gene_dists.1)$distance, as.data.frame(between_gene_dists.2)$distance)
  gene_dists = data.frame(organism = org, distance_to_neighbour = gene_dists)
  all_gene_dists = rbind(all_gene_dists, gene_dists)

  
  ## DISTANCES FROM ASSEMBLED GENES TO ANTISENSE REFERENCE CDS
  org_bs_annot = trinity_translated_bs %>% filter(organism == org_sci)
  org_bs_annot_neg = org_bs_annot %>% 
    filter(base_class == "complete") %>% 
    filter(match_strand == "-")
  
  if(org == "arabidopsis"){org_bs_annot_neg$query_gene = strv_split(org_bs_annot_neg$query_gene, "[.]", 1)}
  
  # find dist to nearest gene on opposite strand
  # split into +/-ve strand genes
  m1 = match(org_bs_annot_neg$query_gene, genes$gene_id)
  # check 'query_gene' matches gene_id, otherwise use 'query'
  if(all(is.na(m1))){
    m2 = match(strv_split(org_bs_annot_neg$query, "[.]", 1), strv_split(ensembl_gff$transcript_id, "[.]", 1))
    org_bs_annot_neg$query_gene = ensembl_gff$gene_id[m2]
  }
  org_bs_annot_neg = org_bs_annot_neg[!is.na(match(org_bs_annot_neg$query_gene, genes$gene_id)),]
  org_bs_annot_neg$gene_strand = strand(genes[match(org_bs_annot_neg$query_gene, genes$gene_id)]) %>% as.character()
  org_bs_annot_neg.posstrand = org_bs_annot_neg[org_bs_annot_neg$gene_strand == "+",]
  org_bs_annot_neg.negstrand = org_bs_annot_neg[org_bs_annot_neg$gene_strand == "-",]
  
  genes.pos = genes[match(org_bs_annot_neg.posstrand$query_gene, genes$gene_id)]
  genes.neg = genes[strand(genes) == "-"]
  pos_dist = GenomicRanges::distanceToNearest(genes.pos, genes.neg, ignore.strand = T)
  
  genes.neg = genes[match(org_bs_annot_neg.negstrand$query_gene, genes$gene_id)]
  genes.pos = genes[strand(genes) == "+"]
  neg_dist = GenomicRanges::distanceToNearest(genes.neg, genes.pos, ignore.strand = T)
  
  org_bs_annot_neg.posstrand$distance_to_neighbour = as.data.frame(pos_dist)$distance
  org_bs_annot_neg.negstrand$distance_to_neighbour = as.data.frame(neg_dist)$distance
  
  # add in distance values (will have NAs for no match)
  org_bs_annot_neg$distance_to_neighbour = org_bs_annot_neg.negstrand$distance_to_neighbour[match(org_bs_annot_neg$subject, org_bs_annot_neg.negstrand$subject)]
  # replace NAs
  org_bs_annot_neg$distance_to_neighbour[org_bs_annot_neg$gene_strand == "+"] = org_bs_annot_neg.posstrand$distance_to_neighbour[match(org_bs_annot_neg$subject[org_bs_annot_neg$gene_strand == "+"], org_bs_annot_neg.posstrand$subject)]
  
  #ggplot(org_bs_annot_neg, aes(x=distance_to_neighbour+0.1, fill=match_type)) + geom_histogram(position="dodge") + scale_x_log10()
  
  negative_annots_with_dists = rbind(negative_annots_with_dists, org_bs_annot_neg)
  
  ## REPLACE -ve complete assembled txs with their +ve version
  replace_with = trinity_translated[trinity_translated$organism == org_sci,]
  
  m = match(org_bs_annot$subject[org_bs_annot$match_strand == "-"], replace_with$subject)
  
  r = m[which(!is.na(m))]
  replaced = plyr::rbind.fill(org_bs_annot[!(org_bs_annot$subject %in% replace_with$subject[r]),], replace_with[r,])
  
  table(replaced$base_class,replaced$match_strand)
  table(org_bs_annot$base_class, org_bs_annot$match_strand)
  
  replaced$organism = paste0(org_sci, " - neg replaced")
  replaced$base_class_f2 = factor(gsub("_", " ", replaced$base_class), levels = c(levels = gsub("_"," ",rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage")))))
  replaced$match_strand = factor(replaced$match_strand, levels=c("+","-"))
  
  trinity_translated_bs = trinity_translated_bs[trinity_translated_bs$organism != paste0(org_sci, " - neg replaced"),]
  trinity_translated_bs = plyr::rbind.fill(trinity_translated_bs, replaced)
  
}

tt_bs_org = strv_split(trinity_translated_bs$organism, "[ ][-]", 1)
tt_org = strv_split(trinity_translated$organism, "[ ][-]", 1)

# make query_gene actually refer to a gene_id, not transcript_id
for(org in c("yeast", "human", "zebrafish", "arabidopsis")){
  org_sci = org_scis[which(orgs==org)]
  
  if(org == "arabidopsis"){
    gff_file = list.files(paste0("data/reference_annotations/", org, "/"), full.names = TRUE)
    gff_file = gff_file[tolower(strv_split(gff_file, "[.]", -1)) %in% c("gff", "gff3", "gtf")]
    gff_file = gff_file[!grepl("utr", gff_file)]
    gtf_type = str_sub(toupper(strv_split(gff_file, "[.]", -1)), 1,3)
    gff_file = gff_file[grep("Araport", gff_file)][1]
    
    ensembl_gff = fread(paste0("grep -v '#' ", gff_file), data.table = F, fill=T, sep='\t')
    colnames(ensembl_gff) = c('seqid','source','type','start','end','score','strand','phase','attributes')
    
    ## IF ARABIDOPSIS
    #ensembl_gff = ensembl_gff[!ensembl_gff$type %in% c("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "protein", "gene","transposable_element",
    #                                                   "transposable_element_gene","transposon_fragment", "pseudogene","pseudogenic_exon", "uORF", "miRNA"),]
    ensembl_gff = ensembl_gff[ensembl_gff$type %in% c("exon"),]
    ensembl_gff$gene_id = strv_split2(ensembl_gff$attributes,"ID[=]", ":")
    
    ensembl_gff$transcript_id = strv_split2(ensembl_gff$attributes,"Parent[=]", "[;]")
    ensembl_gff = separate_rows(ensembl_gff, transcript_id)
    gene_to_tx = ensembl_gff[!duplicated(paste0(ensembl_gff$gene_id, " tx " , ensembl_gff$transcript_id)),]
    
  }else{
    # read in ensembl annotated GFF file
    fa_file = list.files(paste0("data/reference_annotations/", org, "/"), full.names = TRUE)
    fa_file = fa_file[tolower(strv_split(fa_file, "[.]", -1)) %in% c("fa")]
    fa_file  = fa_file[grepl("cdna", fa_file)]
    fa = halpme::read_fasta2df(fa_file)
    fa$seq=NULL
    fa$transcript_id = strv_split(fa$seq_id, "[ ]", 1)
    fa$gene_id = strv_split2(fa$seq_id, "gene[:]", "[ ]")
    
    gene_to_tx = fa
  }
  
  m = match(strv_split(trinity_translated_bs$query[tt_bs_org == org_sci], "[.]", 1), strv_split(gene_to_tx$transcript_id, "[.]", 1))
  trinity_translated_bs$query_gene[tt_bs_org == org_sci] = gene_to_tx$gene_id[m]
  m = match(strv_split(trinity_translated$query[tt_org == org_sci], "[.]", 1), strv_split(gene_to_tx$transcript_id, "[.]", 1))
  trinity_translated$query_gene[tt_org == org_sci] = gene_to_tx$gene_id[m]
}

### WRITE YEAST TT_REPLACED
write.table(trinity_translated_bs[trinity_translated_bs$organism == "s. cerevisiae - neg replaced",], "data/yeast/trinity_translated_negreplaced.txt")


t = table(trinity_translated_bs$organism[which(trinity_translated_bs$base_class_f2 == "complete" & trinity_translated_bs$organism %in% c(org_scis, paste0(org_scis, " - neg replaced")))], 
          trinity_translated_bs$match_strand[which(trinity_translated_bs$base_class_f2 == "complete" & trinity_translated_bs$organism %in% c(org_scis, paste0(org_scis, " - neg replaced")))])

prop_pos = as.data.frame(t / rowSums(t))

prop_pos %>%
  filter(Var2=="+") %>%  
  ggplot(aes(x=Var1, y=Freq)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous("Proportion + strand Blast match (complete ORF match)") + scale_x_discrete("Dataset")

strand_supp_1A = 
  trinity_translated_bs %>% 
  filter(organism %in% c("h. sapiens", "h. sapiens - protein_coding only", "h. sapiens - CDS only",
                         "s. cerevisiae", "s. cerevisiae - 2 samples (116M reads)", "s. cerevisiae - 3 samples (171M reads)")) %>% 
  mutate(org_ordered = factor(gsub("[ ][-][ ]","\n",organism), levels = gsub("[ ][-][ ]","\n",c("h. sapiens", "h. sapiens - protein_coding only", "h. sapiens - CDS only",
                                                                                                "s. cerevisiae", "s. cerevisiae - 2 samples (116M reads)", "s. cerevisiae - 3 samples (171M reads)")))) %>% 
  ggplot(aes(x=match_strand, fill=base_class_f2)) + 
  geom_bar() + facet_wrap(~org_ordered, scales="free", ncol=3, strip.position="right") +
  theme_figure +
  scale_x_discrete("BlastN match strand") + 
  scale_y_continuous("Number of Trinity transcripts") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow=4)) +
  scale_fill_manual(values=base_class_pal2, "Assembled transcript class")


t2 = trinity_translated_bs %>% 
  filter(organism %in% c("h. sapiens", "h. sapiens - protein_coding only", "h. sapiens - CDS only",
                         "s. cerevisiae", "s. cerevisiae - 2 samples (116M reads)", "s. cerevisiae - 3 samples (171M reads)")) %>% 
  mutate(org_ordered = factor(organism, levels = c("h. sapiens", "h. sapiens - protein_coding only", "h. sapiens - CDS only",
                                                   "s. cerevisiae", "s. cerevisiae - 2 samples (116M reads)", "s. cerevisiae - 3 samples (171M reads)")))
t = table(t2$organism[which(t2$base_class_f2 == "complete")], 
          t2$match_strand[which(t2$base_class_f2 == "complete")])

prop_pos = as.data.frame(t / rowSums(t))

strand_supp_1B = prop_pos %>%
  filter(Var2=="+") %>%  
  ggplot(aes(x=Var1, y=Freq)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous("Proportion + strand Blast match\n(complete ORF match)") + scale_x_discrete("Dataset") + 
  theme_figure + theme(panel.grid.major = element_line(), axis.text.x = element_text(angle = 90, hjust = 1))

library(cowplot)
pdf("plots/Supp_Figure_altcds_altdepth.pdf", height=4,width=6.9)
ggdraw() +
  draw_plot(strand_supp_1A, 0,0,0.7,1) + 
  draw_plot(strand_supp_1B, 0.7,0.1,0.3,0.7) +
  draw_plot_label(c("A","B"), size=12, x=c(0,0.7), y=c(1,1))   
dev.off()


strand_supp_3A = 
  trinity_translated %>% 
  mutate(organism = paste0(organism, " - positive only")) %>% 
  plyr::rbind.fill(trinity_translated_bs) %>% 
  filter(organism %in% c(org_scis, paste0(org_scis, " - neg replaced"), paste0(org_scis, " - positive only"))) %>% 
  mutate(org_ordered = gsub("[ ][-][ ]","\n", gsub(" neg ", " negative ", organism))) %>% 
  ggplot(aes(x=match_strand, fill=base_class_f2)) + 
  geom_bar() + facet_wrap(~org_ordered, scales="free_y", ncol=3, strip.position="right") +
  theme_figure +
  scale_x_discrete("BlastN match strand") + 
  scale_y_continuous("Number of Trinity transcripts") +
  #theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol=1)) +
  scale_fill_manual(values=base_class_pal2, "Assembled transcript class") 


t3 = trinity_translated_bs %>% 
  filter(organism %in% c(org_scis, paste0(org_scis, " - neg replaced"), paste0(org_scis, " - positive only"))) %>% 
  mutate(org_ordered = gsub(" neg ", " negative ", organism))

t = table(t3$organism[which(t3$base_class_f2 == "complete")], 
          t3$match_strand[which(t3$base_class_f2 == "complete")])

prop_pos = as.data.frame(t / rowSums(t))

strand_supp_3B = 
  prop_pos %>%
  filter(Var2=="+") %>%  
  ggplot(aes(x=Var1, y=Freq)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous("Proportion + strand Blast match\n(complete ORF match)") + scale_x_discrete("Dataset") + 
  theme_figure + theme(panel.grid.major = element_line(), axis.text.x = element_text(angle = 90, hjust = 1))

leg <- ggplot_legend(strand_supp_3A)
pdf("plots/Supp_Figure_negativereplacement.pdf", height=6,width=6.9)
ggdraw() +
  draw_plot(strand_supp_3A + theme(legend.position = "none"), 0,0,0.7,1) + 
  draw_plot(strand_supp_3B, 0.7,0,0.3,0.65) +
  draw_plot(leg, 0.7,0.7,0.3,0.3) +
  draw_plot_label(c("A","B"), size=12, x=c(0,0.7), y=c(1,0.65))   
dev.off()



all_gene_dists$organism = org_scis[match(all_gene_dists$organism, orgs)]
dist_summary = aggregate(distance_to_neighbour ~ organism, all_gene_dists, median)
prop_overlap = table(all_gene_dists$organism,all_gene_dists$distance_to_neighbour == 0)
dist_summary = as.data.frame(prop_overlap / rowSums(prop_overlap)) %>% dplyr::rename(organism=Var1,prop_overlapping=Freq) %>% 
  filter(Var2 == TRUE) %>% left_join(dist_summary,., by=c('organism'))

dist_summary$y_val = 2700
dist_summary$y_val[dist_summary$organism == "h. sapiens"] = 2200

dist_summary = arrange(dist_summary, desc(prop_overlapping))
dist_summary$y_val_annot = c(1800,1600,1400,1200)
dist_summary$annot = paste0(dist_summary$organism, ": ", round(dist_summary$prop_overlapping, 3)*100, "%")


reference_dist_p1 = ggplot(all_gene_dists, 
                           aes(x=distance_to_neighbour+0.1, fill=organism)) + scale_y_continuous("Count")+ 
  geom_histogram(position='dodge', bins=50) + scale_x_log10("Distance (nt) to nearest opposite strand gene + 0.1")  + 
  geom_vline(data=dist_summary, aes(xintercept = distance_to_neighbour, col=organism), linetype=2) + 
  #annotate(data=dist_summary, geom="text", aes(x=distance_to_neighbour, label=distance_to_neighbour)) +
  geom_label(data=dist_summary,aes(x=distance_to_neighbour, label=distance_to_neighbour, y = y_val, fill=organism), alpha=0.2, hjust=0, size=2) + 
  geom_text(data=dist_summary,aes(x=0.5+0.1, label=paste0(organism, ": ", round(prop_overlapping, 3)*100 , "%"), y = y_val_annot, col=organism), hjust=0, size=2) + 
  annotate("text", x=0.5+0.1, y= 2200, hjust=0, label="Percent overlapping\nCDS on opposite strand", size=2)+
  annotate("rect", xmin=0+0.08, xmax=0.2+0.1, ymin=0, ymax=2000, alpha=0.1) +
  #annotate("segment", x=0.2+0.1, y=2000, xend=0.5+01, yend=2200) + 
  theme_figure + scale_fill_brewer(palette = "Set2", "Species") + scale_color_brewer(palette = "Set2", "Species")


dist_summary_real = aggregate(distance_to_neighbour ~ organism, negative_annots_with_dists, median)
prop_overlap = table(negative_annots_with_dists$organism,negative_annots_with_dists$distance_to_neighbour == 0)
dist_summary_real = as.data.frame(prop_overlap / rowSums(prop_overlap)) %>% dplyr::rename(organism=Var1, prop_overlapping=Freq) %>% 
  filter(Var2 == TRUE) %>% left_join(dist_summary_real,., by=c('organism'))

dist_summary_real = as.data.frame(prop_overlap) %>% dplyr::rename(organism=Var1, prop_overlapping=Freq) %>% 
  filter(Var2 == TRUE) %>% left_join(dist_summary_real,., by=c('organism'), suffix=c("",".real"))

negative_annots_with_dists$match_type2 = negative_annots_with_dists$match_type
negative_annots_with_dists$match_type2 = str_replace(negative_annots_with_dists$match_type2, "multi_hit", "Multiple hits") %>% str_replace("single_hit", "Single hit")

#max_heights = aggregate(log(distance_to_neighbour+0.1)~organism, negative_annots_with_dists, function(x)  max(table(cut(x, breaks=30))))
max_heights_func = function(distance_vector, max_zlog=max){
  z_log = log10(distance_vector + 0.1)
  z_log = sort(z_log)
  z_log_widths = (max_zlog-(min(z_log))) / 30
  bins = (min(z_log)) + ((z_log_widths)*(0:30))
  heights = vector()
  for(i in 1:(length(bins)-1)){
    heights[i] = length(which(z_log >= bins[i] & z_log < bins[i+1]))
  }
  return(max(heights))
}



max_heights = aggregate(distance_to_neighbour~organism, negative_annots_with_dists, 
                        function(x) max_heights_func(x,max_zlog=max(log10(negative_annots_with_dists$distance_to_neighbour+0.1))))
dist_summary_real$y_max = max_heights[match(dist_summary_real$organism, max_heights$organism),2]

real_dists_p2 = ggplot(negative_annots_with_dists) + 
  geom_histogram(position="stack", aes(x=distance_to_neighbour+0.1, fill=match_type2)) + facet_wrap(~organism, scales="free_y", ncol=4) + 
  geom_vline(data=dist_summary_real, aes(xintercept = distance_to_neighbour), linetype=2) + 
  geom_text(data=dist_summary_real, aes(x=0.2+0.1, label=paste0(round(prop_overlapping, 3)*100, "%"), y=prop_overlapping.real), hjust=0, vjust=1, size=2) + 
  geom_label(data=dist_summary_real,aes(x=distance_to_neighbour, label=distance_to_neighbour, y = y_max), alpha=0.2, hjust=0, size=2) + 
  #ggpmisc::geom_label_npc(data=dist_summary_real,aes(label=distance_to_neighbour), alpha=0.2, hjust=1, npcx=0.5, npcy=1) +
  scale_x_log10("Distance (nt) to nearest opposite strand gene + 0.1") + theme_figure + scale_y_continuous("Count") + 
  scale_fill_discrete("BlastN\nmatch type")


pdf("plots/Supp_Figure_genedists.pdf", height=4,width=6.9)
ggdraw() +
  draw_plot(reference_dist_p1 + ggtitle("Reference"), 0,0.5,1,0.5) + 
  draw_plot(real_dists_p2 + ggtitle("Assembled"), 0,0,1,0.5) +
  draw_plot_label(c("A","B"), size=12, x=c(0,0), y=c(1,0.5))   
dev.off()


tt_neg = 
trinity_translated_bs %>% 
  filter(organism %in% c(paste0(org_scis, " - neg replaced"))) %>% 
  mutate(org_ordered = gsub("[ ][-][ ]","\n", gsub(" neg ", " negative ", organism)))
 
td_all = NULL
tt_all = NULL
for(org in c("human", "arabidopsis", 'zebrafish', 'yeast')){
  
  tt = tt_neg[tt_neg$org==org,]
  load(paste0("data/", org, "/processed_ensembl.Rdata"))
  
  ensembl_peps$gene_id = NA
  ensembl_peps$gene_id[grepl("gene[:]", ensembl_peps$seq_id)] = strv_split2(ensembl_peps$seq_id[grepl("gene[:]", ensembl_peps$seq_id)], "gene[:]", "[ ]")
  if(org == "arabidopsis"){
    ensembl_peps$gene_id = ensembl_gff_cds$gene_id[match(ensembl_peps$transcript, ensembl_gff_cds$transcript_id)]
  }
  borf_results = fread(paste0("data/", org, "/borf/Trinity.txt"), data.table = F)
  borf_seq = halpme::read_fasta2df(paste0("data/", org, "/borf/Trinity.pep"))
  
  borf_seq = borf_seq[borf_results$transcript_id %in% tt$subject,]
  borf_results = borf_results[borf_results$transcript_id %in% tt$subject,]
  borf_results$borf_aa = borf_seq$seq[match(borf_results$orf_id, paste0(">", borf_seq$seq_id))]
  
  borf_blast = fread(paste0("data/", org, "/borf/borf_blastp.txt"), data.table = F)
  colnames(borf_blast) = c("query", "subject", "p_ident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","sframe", "qframe")
  borf_blast$subject_transcript_id = ensembl_peps$transcript[match(borf_blast$subject, strv_split(ensembl_peps$seq_id, "[ ]", 1))]
  borf_blast$subject_gene_id = ensembl_peps$gene_id[match(borf_blast$subject, strv_split(ensembl_peps$seq_id, "[ ]", 1))]
  
  borf_results$base_hit = tt$query[match(borf_results$transcript_id, tt$subject)]
  borf_results$hit_blastp = borf_blast$subject_transcript_id[match(gsub(">","",strv_split(borf_results$orf_id, "[ ]",1)), borf_blast$query)]
  borf_results$base_in_blastp = ifelse(paste0(gsub(">","",strv_split(borf_results$orf_id, "[ ]",1)), borf_results$base_hit) %in% paste0(borf_blast$query, borf_blast$subject_transcript_id), 1, 0)
  borf_results$base_hit_gene = ensembl_peps$gene_id[match(borf_results$base_hit, ensembl_peps$transcript)]
  
  borf_results$base_in_blastp_gene = ifelse(paste0(gsub(">","",strv_split(borf_results$orf_id, "[ ]",1)), borf_results$base_hit_gene) %in% paste0(borf_blast$query, borf_blast$subject_gene_id), 1, 0)
  
  
  colnames(borf_blast) = paste0("blastp_borf_",colnames(borf_blast))
  borf_results = cbind(borf_results, borf_blast[match(gsub(">","",strv_split(borf_results$orf_id, "[ ]",1)), borf_blast$blastp_borf_query),])
  
  # write negative/positive strand TDs for interproscan
  seqs_b = borf_results[, c(1,2)]
  seqs_b$seq_id = paste0(">", seqs_b$seq_id)
  seqs_b$seq = gsub("[*]","", seqs_b$seq)
  seqs_b_pos = seqs_b[borf_results$strand =="+",]

  transdecoder_results = halpme::read_fasta2df(paste0("data/", org, "/Trinity.fasta.transdecoder_dir/longest_orfs.pep"))
  transdecoder_results$orf_id = strv_split(transdecoder_results$seq_id, "[ ]", 1)
  transdecoder_results$transcript_id = strv_split(transdecoder_results$orf_id, "[.]p", 1)
  transdecoder_results$strand = strv_split2(transdecoder_results$seq_id, "[(]", "[)]")
  
  locs = strv_split(transdecoder_results$seq_id, "[ ]", 5)
  transdecoder_results$start_site_nt = as.numeric(strv_split2(locs, "[:]", "[-]"))
  transdecoder_results$stop_site_nt = as.numeric(strv_split2(locs, "[-]", "[(]"))
  transdecoder_results$orf_length_aa = as.numeric(strv_split2(transdecoder_results$seq_id, "len[:]", "[ ]"))
  transdecoder_results$orf_length_aa = nchar(gsub("[*]", "", transdecoder_results$seq))
  transdecoder_results$orf_class = strv_split2(transdecoder_results$seq_id, "type[:]", "[ ]")
  
  transdecoder_results = transdecoder_results[transdecoder_results$transcript_id %in% tt$subject,]
  transdecoder_results = arrange(transdecoder_results, desc(orf_length_aa))
  transdecoder_results$org =org
  
  
  td_blast = fread(paste0("data/", org, "/Trinity.fasta.transdecoder_dir/transdecoder_blastp.txt"), data.table = F)
  
  
  colnames(td_blast) = c("query", "subject", "p_ident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","sframe", "qframe")
  td_blast$subject_transcript_id = ensembl_peps$transcript[match(td_blast$subject, strv_split(ensembl_peps$seq_id, "[ ]", 1))]
  td_blast$subject_gene_id = ensembl_peps$gene_id[match(td_blast$subject, strv_split(ensembl_peps$seq_id, "[ ]", 1))]
  
  transdecoder_results$base_hit = tt$query[match(transdecoder_results$transcript_id, tt$subject)]
  transdecoder_results$hit_blastp = td_blast$subject_transcript_id[match(transdecoder_results$orf_id, td_blast$query)]
  transdecoder_results$base_in_blastp = ifelse(paste0(transdecoder_results$orf_id, transdecoder_results$base_hit) %in% paste0(td_blast$query, td_blast$subject_transcript_id), 1, 0)
  
  transdecoder_results$base_hit_gene = ensembl_peps$gene_id[match(transdecoder_results$base_hit, ensembl_peps$transcript)]
  transdecoder_results$base_in_blastp_gene = ifelse(paste0(transdecoder_results$orf_id, transdecoder_results$base_hit_gene) %in% paste0(td_blast$query, td_blast$subject_gene_id), 1, 0)
  
  
  colnames(td_blast) = paste0("blastp_td_",colnames(td_blast))
  transdecoder_results = cbind(transdecoder_results, td_blast[match(transdecoder_results$orf_id, td_blast$blastp_td_query),])
  
  transdecoder_results_add = cbind(transdecoder_results, tt[match(transdecoder_results$transcript_id, tt$subject),])
  
  td_all = rbind(td_all, transdecoder_results_add)
  transdecoder_results$org =NULL
  
  transdecoder_results = transdecoder_results[transdecoder_results$orf_length_aa >=100,]
  #transdecoder_results = transdecoder_results[transdecoder_results$strand == "+",]
  transdecoder_results = arrange(transdecoder_results, desc(orf_length_aa))
  transdecoder_results = transdecoder_results[!duplicated(transdecoder_results$transcript_id),]
  
  # write negative/positive strand TDs for interproscan
  seqs_td = transdecoder_results[, c(1,2)]
  seqs_td$seq_id = paste0(">", seqs_td$seq_id)
  seqs_td$seq = gsub("[*]","", seqs_td$seq)
  seqs_td_pos = seqs_td[transdecoder_results$strand =="+",]
  seqs_td_neg = seqs_td[transdecoder_results$strand =="-",]
  
  colnames(borf_results) = paste0(colnames(borf_results), "_borf")
  #tt = left_join(tt, borf_results, by=c('subject'='transcript_id_borf'))
  tt = cbind(tt, borf_results[match(tt$subject, borf_results$transcript_id_borf),])
  
  
  colnames(transdecoder_results) = paste0(colnames(transdecoder_results), "_transdecoder")
  tt = cbind(tt, transdecoder_results[match(tt$subject, transdecoder_results$transcript_id_transdecoder),])
  
  tt$orf_class_transdecoder[is.na(tt$orf_class_transdecoder)] = "no_ORF"
  tt$orf_class_transdecoder[which(tt$orf_class_transdecoder == "3prime_partial")] = "incomplete_3prime"
  tt$orf_class_transdecoder[which(tt$orf_class_transdecoder == "5prime_partial")] = "incomplete_5prime"
  tt$orf_class_transdecoder[which(tt$orf_class_transdecoder == "internal")] = "incomplete"
  
  tt$orf_class_borf[is.na(tt$orf_class_borf)] = "no_ORF"
  tt$org=org
  tt_all=plyr::rbind.fill(tt_all, tt)
  
}

tt_all$same_as_base = ifelse(tt_all$base_class == tt_all$orf_class_transdecoder & tt_all$base_class == tt_all$orf_class_borf & tt_all$base_in_blastp_borf == 1 & tt_all$base_in_blastp_transdecoder== 1, "same", "no")
tt_all$same_as_base[which((tt_all$base_class == tt_all$orf_class_transdecoder & tt_all$base_in_blastp_transdecoder== 1) & (tt_all$base_class != tt_all$orf_class_borf | tt_all$base_in_blastp_borf != 1))] = "td_correct"
tt_all$same_as_base[which((tt_all$base_class != tt_all$orf_class_transdecoder | tt_all$base_in_blastp_transdecoder!= 1) & (tt_all$base_class == tt_all$orf_class_borf & tt_all$base_in_blastp_borf == 1))] = "borf_correct"
tt_all$same_as_base[which((tt_all$base_class != tt_all$orf_class_transdecoder | tt_all$base_in_blastp_transdecoder!= 1) & (tt_all$base_class != tt_all$orf_class_borf | tt_all$base_in_blastp_borf != 1))] = "both_incorrect"
tt_all$same_as_base = factor(tt_all$same_as_base, levels = rev(c("same", "td_correct", "borf_correct", "both_incorrect")))
tt_all$same_as_base_c = as.character(tt_all$same_as_base)
tt_all$same_as_base_c[tt_all$same_as_base_c == "same"] = "Borf and Transdecoder"
tt_all$same_as_base_c[tt_all$same_as_base_c == "td_correct"] = "Transdecoder"
tt_all$same_as_base_c[tt_all$same_as_base_c == "borf_correct"] = "Borf"
tt_all$same_as_base_c[tt_all$same_as_base_c == "both_incorrect"] = "Neither"
tt_all$same_as_base_c = factor(tt_all$same_as_base_c, levels = rev(c("Borf and Transdecoder", "Transdecoder", "Borf", "Neither")))


# incorrect vs. correct

keephigh = which(tt_all$base_class %in% c("complete", "incomplete", "incomplete_5prime","incomplete_3prime","insufficient_coverage") & 
                   tt_all$lcs_len >=100 & tt_all$read_counts >= 100)

keephigh10 = which(tt_all$base_class %in% c("complete", "incomplete", "incomplete_5prime","incomplete_3prime","insufficient_coverage") & 
                     tt_all$lcs_len >=100 & tt_all$read_counts >= 10)
plot_corrincorr = 
  ggplot(tt_all[keephigh,], aes(x=base_class, fill=same_as_base_c)) + geom_bar() + 
  facet_wrap(~organism, scales="free_x", ncol=1) + 
  scale_fill_manual(values=c("#fb9a99", "#00B8D9","#FFAB00","#b2df8a"), "Predicted ORF blasts\nto correct gene")  + theme_figure + 
  scale_x_discrete("Standard ORF class") + scale_y_continuous("Transcripts") + coord_flip() + theme(legend.position = c(0.8,0.95)) + 
  guides(fill=guide_legend(ncol=2)) + ggtitle("Predicted ORF type comparison to base type")

table(tt_all_exp$base_in_blastp_borf[tt_all_exp$base_class == "insufficient_coverage"], 
      tt_all_exp$org[tt_all_exp$base_class == "insufficient_coverage"])
table(tt_all_exp$base_in_blastp_transdecoder[tt_all_exp$base_class == "insufficient_coverage"], 
      tt_all_exp$org[tt_all_exp$base_class == "insufficient_coverage"])


# method accuracy

totals = table(tt_all$same_as_base[keephigh], tt_all$organism[keephigh]) %>% as.data.frame()
total_ts = aggregate(Freq ~ Var2, totals, sum)
total_ts$borf_correct = aggregate(Freq ~ Var2, totals[totals$Var1 %in% c("same", "borf_correct"),], sum)[,2]
total_ts$td_correct = aggregate(Freq ~ Var2, totals[totals$Var1 %in% c("same", "td_correct"),], sum)[,2]
total_ts$borf_acc = total_ts$borf_correct / total_ts$Freq
total_ts$td_acc = total_ts$td_correct / total_ts$Freq

total_ts_long = data.frame(organism = rep(total_ts$Var2, 2), 
                           total_transcripts = rep(total_ts$Freq, 2),
                           correct_n = c(total_ts$borf_correct, total_ts$td_correct),
                           accuracy = c(total_ts$borf_acc, total_ts$td_acc),
                           method = rep(c("Borf", "Transdecoder"), each=4))

plot_method_accuracy = ggplot(total_ts_long, aes(x=organism, y=accuracy, col=method)) + 
  geom_point() + 
  theme_figure + 
  scale_color_manual(values = c("#00B8D9", "#FFAB00"), "Method") + 
  ggtitle("Prediction accuracy") + 
  scale_x_discrete("Species") + scale_y_continuous("Accuracy") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

colnames(total_ts)[1:2] = c("Species/Annotation", "Total transcripts")
write_csv(total_ts, "supp_negreplace_accuracy.csv")

pdf("plots/Supp_Figure_negativereplacementaccuracy.pdf", height=3, width=3.3)
plot_method_accuracy
dev.off()

