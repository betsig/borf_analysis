### analyse / plot trinity butts

library(data.table)
library(tidyverse)
library(stringr)
library(halpme)
library(GenomicRanges)
library(stringdist)
library(Biostrings)
devtools::install_github("jonocarroll/ggeasy")
library(ggeasy)



options(stringsAsFactors = F)

at_color_pals = read.delim("at_colour_pal.txt")

base_class_pal_muted = c("#fb8072","#fccde5","#ffffb3","#d9d9d9","#ccebc5","#bc80bd","#ffed6f","#80b1d3","#b3de69")

base_class_pal = c("#e31a1c","#cab2d6","#fdbf6f","#a6cee3","#b2df8a","#6a3d9a","#ff7f00","#1f78b4","#33a02c")

species_pal = c("#36B37E","#00B8D9","#FF5630","#FFAB00")

inten = 300
base_class_pal2 = c(at_color_pals$hex_value[at_color_pals$colour == "N" & at_color_pals$intensity ==inten-250],
                    at_color_pals$hex_value[at_color_pals$colour == "N" & at_color_pals$intensity ==inten+100],
                    at_color_pals$hex_value[at_color_pals$colour == "P" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "B" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "T" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "G" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "Y" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "R" & at_color_pals$intensity ==inten])

base_class_pal2.1 = c(at_color_pals$hex_value[at_color_pals$colour == "N" & at_color_pals$intensity ==inten-250],
                      at_color_pals$hex_value[at_color_pals$colour == "N" & at_color_pals$intensity ==inten-100],
                    at_color_pals$hex_value[at_color_pals$colour == "N" & at_color_pals$intensity ==inten+200],
                    at_color_pals$hex_value[at_color_pals$colour == "P" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "G" & at_color_pals$intensity ==inten-100],
                    at_color_pals$hex_value[at_color_pals$colour == "B" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "T" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "G" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "Y" & at_color_pals$intensity ==inten],
                    at_color_pals$hex_value[at_color_pals$colour == "R" & at_color_pals$intensity ==inten])


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

ggplot_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

trinity_translated = NULL
for(org in c('yeast', 'arabidopsis', 'zebrafish','human')){
    ## BOTH STRANDS DATA
    tt = fread(paste0('data/', org, "/trinity_translated_both_strands.txt"), data.table = F)
    trinity_seqs = fread(paste0("data/", org, "/trinity_with_kallisto.txt"), data.table = F)
    colnames(trinity_seqs) = c("full_id", "subject", "trinity_gene", "trinity_length_nt", "read_counts", "reads_per_kb")
    trinity_seqs$base_class = "no_hits"
    tt = plyr::rbind.fill(tt, trinity_seqs[which(!(trinity_seqs$subject %in% tt$subject)),-1])
    tt$org = org
    
    trinity_translated = plyr::rbind.fill(trinity_translated, tt)
    
}


trinity_translated$base_class[trinity_translated$read_counts < 10] = "<10 reads"
trinity_translated$base_class_factor = factor(trinity_translated$base_class, levels = rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage", "no_hits", "<10 reads")))

trinity_translated$organism = "h. sapiens"
trinity_translated$organism[trinity_translated$org == "zebrafish"] = "d. rerio"
trinity_translated$organism[trinity_translated$org == "arabidopsis"] = "a. thaliana"
trinity_translated$organism[trinity_translated$org == "yeast"] = "s. cerevisiae"

trinity_translated$base_class_f2 = factor(gsub("_", " ", trinity_translated$base_class), levels = c(levels = gsub("_"," ",rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage", "no_hits", "<10 reads")))))
table(trinity_translated$org)


# make annotation df for displaying percentages on plot
tt_class_totals = table(trinity_translated$base_class_f2 , trinity_translated$organism) %>% as.data.frame() %>% filter(Freq>0)
tt_totals = table(trinity_translated$organism) %>% as.data.frame() %>% filter(Freq>0)
tt_class_totals = tt_class_totals %>% left_join(tt_totals, by=c('Var2'='Var1')) %>% mutate(percent =paste0(round((Freq.x/Freq.y)*100), "%"))
colnames(tt_class_totals)[1:2] = c("base_class_f2", "organism")
tt_class_totals$y_val = 0
tt_class_totals$y_val[tt_class_totals$base_class_f2=="complete"] = tt_class_totals$Freq.y[tt_class_totals$base_class_f2=="complete"] * 0.1
tt_class_totals$y_val[tt_class_totals$base_class_f2=="complete" & tt_class_totals$organism == "h. sapiens"] = 
  tt_class_totals$Freq.y[tt_class_totals$base_class_f2=="complete" & tt_class_totals$organism == "h. sapiens"] * 0.028

tt_class_totals$y_val[tt_class_totals$base_class_f2=="<10 reads"] = tt_class_totals$Freq.y[tt_class_totals$base_class_f2=="<10 reads"] * 0.9
new_tots = aggregate(Freq.x ~ organism, tt_class_totals[!(tt_class_totals$base_class_f2 %in% c("no hits", "<10 reads")),], sum)
m = match(tt_class_totals$organism, new_tots$organism)
tt_class_totals$y_val[tt_class_totals$base_class_f2=="no hits"] = (new_tots$Freq.x[m] + (tt_class_totals$Freq.x*0.5))[tt_class_totals$base_class_f2=="no hits"] 


tt_class_totals$text_col = "none"
tt_class_totals$text_col[tt_class_totals$base_class_f2=="complete"] = "white"
tt_class_totals$text_col[tt_class_totals$base_class_f2=="no hits"] = "white"
tt_class_totals$text_col[tt_class_totals$base_class_f2=="<10 reads"] = "black"

write_delim(tt_class_totals[,c(2,1,3,5)], "data/supp_table_assembled_classes.tsv", delim = "\t")

# FIGURE 1A- number of assembled transcripts in each class - i.e. most are garbage.
plot_n_trans =
  ggplot(trinity_translated, 
       aes(x=organism, fill=base_class_f2)) + 
  geom_bar() + 
  facet_wrap(~organism,ncol=1,  scales='free') + 
  coord_flip() + 
  theme_figure + 
  ggeasy::easy_remove_y_axis() + 
  theme(legend.position = "bottom", axis.text.x =element_text(hjust=1))   + 
  guides(fill = guide_legend(nrow=4, title.position = "top", title.hjust = 0.5)) +
  scale_y_continuous("Number of Trinity transcripts", ) +
  guides(col = guide_legend(nrow = 4)) +
  scale_fill_manual(values=base_class_pal2.1, "Assembled transcript class")+ 
  ggtitle("Trinity transcripts assembled")+
  geom_text(data=tt_class_totals[tt_class_totals$base_class_f2 %in% c("complete", "no hits", "<10 reads"),], aes(label=percent, y=y_val, col=text_col), size=2) + 
  scale_color_manual(values=c("black","white")) + guides(color=FALSE)
####

ret_keep = which(trinity_translated$read_counts >=100 & 
                   trinity_translated$lcs_len >= 10 &
                   trinity_translated$base_class %in% c("complete", "incomplete", "incomplete_5prime","incomplete_3prime"))

tt_ret_class_totals = table(trinity_translated$base_class_f2[ret_keep] , trinity_translated$organism[ret_keep]) %>% as.data.frame() %>% filter(Freq>0)
tt_ret_totals = table(trinity_translated$organism[ret_keep]) %>% as.data.frame() %>% filter(Freq>0)
tt_ret_class_totals = tt_ret_class_totals %>% left_join(tt_ret_totals, by=c('Var2'='Var1')) %>% mutate(percent =paste0(round((Freq.x/Freq.y)*100), "%"))
colnames(tt_ret_class_totals)[1:2] = c("base_class_f2", "organism")
tt_ret_class_totals$y_val = 0
tt_ret_class_totals$y_val[tt_ret_class_totals$base_class_f2=="complete"] = tt_ret_class_totals$Freq.y[tt_ret_class_totals$base_class_f2=="complete"] * 0.1
tt_ret_class_totals$text_col = "none"
tt_ret_class_totals$text_col[tt_ret_class_totals$base_class_f2=="complete"] = "white"
tt_ret_class_totals$text_col[tt_ret_class_totals$base_class_f2=="no hits"] = "black"

# FIGURE 1C 
plot_n_trans_retained=
  ggplot(trinity_translated[which(trinity_translated$read_counts >=100 & 
                                    trinity_translated$lcs_len >= 10 &
                                    trinity_translated$base_class %in% c("complete", "incomplete", "incomplete_5prime","incomplete_3prime")),], 
         aes(x=organism, fill=base_class_f2)) + 
  geom_bar() + 
  facet_wrap(~organism,ncol=1,scales='free') + 
  coord_flip() + 
  theme_figure + 
  ggeasy::easy_remove_y_axis() + 
  easy_move_legend("bottom") + 
  scale_y_continuous("Number of Trinity transcripts") +
  guides(col = guide_legend(nrow = 2, ncol=2)) +
  theme(axis.text.x =element_text(hjust=1))   + 
  scale_fill_manual(values=base_class_pal2[c(5:8)], "Assembled transcript class")+ ggtitle("Trinity transcripts passing filters") +
  geom_text(data=tt_ret_class_totals[tt_ret_class_totals$base_class_f2 %in% c("complete"),], aes(label=percent, y=y_val), size=2, col="white")

## gene_coverage
ens_gene_cov=read.delim("data/ensembl_pcgene_coverage_bothstrands.txt")
ens_gene_cov$class = factor(ens_gene_cov$class, levels = rev(c("assembled_count_100_ORF_match","assembled_count_100", "assembled_count_10_ORF_match","assembled_count_10", "not_assembled")))

ens_gene_cov$class2 = as.character(ens_gene_cov$class)
ens_gene_cov$class2[ens_gene_cov$class == "not_assembled"] = "not assembled"
ens_gene_cov$class2[ens_gene_cov$class == "assembled_count_10"] = "counts >= 10"
ens_gene_cov$class2[ens_gene_cov$class == "assembled_count_10_ORF_match"] = "counts >= 10 & ORF match"
ens_gene_cov$class2[ens_gene_cov$class == "assembled_count_100"] = "counts >= 100"
ens_gene_cov$class2[ens_gene_cov$class == "assembled_count_100_ORF_match"] = "counts >= 100 & ORF match"

ens_gene_cov$class2 = factor(ens_gene_cov$class2, levels = rev(c("counts >= 100 & ORF match","counts >= 100", "counts >= 10 & ORF match","counts >= 10", "not assembled")))


table(ens_gene_cov$organism, ens_gene_cov$class2) / rowSums(table(ens_gene_cov$organism, ens_gene_cov$class2) )


ens_gene_cov_class_totals = table(ens_gene_cov$class2 , ens_gene_cov$organism) %>% as.data.frame() %>% filter(Freq>0)
ens_gene_cov_totals = table(ens_gene_cov$organism) %>% as.data.frame() %>% filter(Freq>0)
ens_gene_cov_class_totals = ens_gene_cov_class_totals %>% left_join(ens_gene_cov_totals, by=c('Var2'='Var1')) %>% mutate(percent =paste0(round((Freq.x/Freq.y)*100), "%"))
colnames(ens_gene_cov_class_totals)[1:2] = c("class2", "organism")
ens_gene_cov_class_totals$y_val = 0
ens_gene_cov_class_totals$y_val[ens_gene_cov_class_totals$class2=="counts >= 100 & ORF match"] = ens_gene_cov_class_totals$Freq.y[ens_gene_cov_class_totals$class2=="counts >= 100 & ORF match"] * 0.1
ens_gene_cov_class_totals$y_val[ens_gene_cov_class_totals$class2=="not assembled"] = ens_gene_cov_class_totals$Freq.y[ens_gene_cov_class_totals$class2=="not assembled"] * 0.9
ens_gene_cov_class_totals$text_col = "none"
ens_gene_cov_class_totals$text_col[ens_gene_cov_class_totals$class2=="counts >= 100 & ORF match"] = "white"
ens_gene_cov_class_totals$text_col[ens_gene_cov_class_totals$class2=="not assembled"] = "black"

plot_ens_gene_cov = 
  ggplot(ens_gene_cov, aes(x=organism, fill=class2)) +
  geom_bar() + 
  facet_wrap(~organism,ncol=1, nrow=4, scales='free') + 
  coord_flip() + 
  theme_figure + 
  scale_y_continuous("Reference protein coding genes") +
  scale_fill_manual("Trinity transcript match",values= rev(at_color_pals$hex_value[at_color_pals$colour == "N"])[seq(6, length(which(at_color_pals$colour == "N")), 3)]) +
  ggeasy::easy_remove_y_axis() + 
  theme(legend.position = "bottom", axis.text.x =element_text(hjust=1))   + 
  guides(fill = guide_legend(nrow=3, title.position = "top", title.hjust = 0.5)) +
  ggtitle("Reference transcriptome coverage")+
  geom_text(data=ens_gene_cov_class_totals[ens_gene_cov_class_totals$class2 %in% c("counts >= 100 & ORF match", "not assembled"),], aes(label=percent, y=y_val, col=text_col), size=2) + 
  scale_color_manual(values=c("black","white")) + guides(color=FALSE)


plot_ens_orf_cov = 
ggplot(trinity_translated[which(trinity_translated$base_class != "insufficient_coverage"& 
                                  trinity_translated$lcs_len >= 10 &
                                  trinity_translated$base_class != "complete" & trinity_translated$read_counts >= 100),], 
       aes(y=lcs_percent_of_ens, x=organism,fill=organism)) + geom_boxplot(outlier.shape = NA, size=0.25)+ 
  scale_fill_brewer(palette = "Set2") + 
  scale_y_continuous("Percent coverage of reference ORF")+
  theme_figure +scale_x_discrete("Species") +
  easy_remove_legend() + ggtitle("Partially assembled\ntranscript\nORF coverage (%)") + theme(axis.text.x=element_text(angle=90,hjust=1))

plot_ens_orf_cov_aa = 
ggplot(trinity_translated[which(trinity_translated$base_class != "insufficient_coverage"& 
                                  trinity_translated$lcs_len >= 10 &
                                  trinity_translated$base_class != "complete" & trinity_translated$read_counts >= 100),], 
       aes(y=lcs_diff_to_ens, x=organism,fill=organism)) + geom_boxplot(outlier.shape = NA, size=0.25)+ 
  scale_fill_brewer(palette = "Set2") + 
  scale_y_log10("AA difference to reference ORF")+
  theme_figure +
  scale_x_discrete("Species") +
  easy_remove_legend() + ggtitle("Partially assembled\ntranscript\nORF coverage (AA)") + theme(axis.text.x=element_text(angle=90,hjust=1))

aggregate(lcs_diff_to_ens ~ org, trinity_translated[which(trinity_translated$base_class != "insufficient_coverage"& 
                                                            trinity_translated$lcs_len >= 10 &
                                                            trinity_translated$base_class != "complete" & trinity_translated$read_counts >= 100),], median )

aggregate(lcs_percent_of_ens ~ org, trinity_translated[which(trinity_translated$base_class != "insufficient_coverage"& 
                                                               trinity_translated$lcs_len >= 10 &
                                                               trinity_translated$base_class != "complete" & trinity_translated$read_counts >= 100),], median )


library(cowplot)
## FIGURE 1
leg_1 = ggplot_legend(plot_n_trans)
leg_2 = ggplot_legend(plot_ens_gene_cov)

pdf("plots/Figure1.pdf", height=6,width=6.6)
ggdraw() +
  draw_plot(plot_n_trans+ theme(legend.position = "none"), 0,0.55,0.5,0.45) + 
  draw_plot(plot_ens_gene_cov + theme(legend.position = "none"), 0.5,0.55,0.5,0.45)+
  draw_plot(leg_1, 0,0.45,0.5,0.1)+
  draw_plot(leg_2, 0.5,0.45,0.5,0.1)+
  
  draw_plot(plot_n_trans_retained + theme(legend.position = "none"), 0,0,0.5,0.45)+
  draw_plot(plot_ens_orf_cov + theme(legend.position = "none"), 0.5,0,0.25,0.45)+
  draw_plot(plot_ens_orf_cov_aa + theme(legend.position = "none"), 0.75,0,0.25,0.45)+
  #draw_plot(leg, 0.86,0,0.14,1) +
  draw_plot_label(c("A","B","C","D","E"), size=12, x=c(0,0.5,0,0.5,0.75), y=c(1,1,0.45,0.45, 0.45))   
dev.off()


trinity_translated$base_class_5 = NA
trinity_translated$base_class_5[trinity_translated$base_class %in% c("complete", "incomplete_3prime")] = "complete"
trinity_translated$base_class_5[trinity_translated$base_class %in% c("incomplete", "incomplete_5prime")] = "incomplete_5prime"

trinity_translated$base_class_5 = factor(trinity_translated$base_class_5, levels = c("incomplete_5prime", "complete"))
trinity_translated$base_class_tag = factor(trinity_translated$base_class_tag, levels = c("partial_match", "full_match"))

com_inc_10 = ggplot(trinity_translated[trinity_translated$read_counts >=10 & trinity_translated$base_class_5 %in% c("complete", "incomplete_5prime"),], 
                    aes(x=base_class_5, fill=base_class_tag)) + geom_bar() + facet_wrap(~organism, ncol=4, scales="free_y")+ 
  scale_fill_manual("Blast match type", values = rev(at_color_pals$hex_value[at_color_pals$colour == "N" & at_color_pals$intensity %in% c(800, 100)])) + 
  theme_figure + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(">=10 reads")  + scale_x_discrete("Class")+ scale_y_continuous("Count")

com_inc_100 = ggplot(trinity_translated[trinity_translated$read_counts >=100 & trinity_translated$base_class_5 %in% c("complete", "incomplete_5prime"),], 
                     aes(x=base_class_5, fill=base_class_tag)) + geom_bar() + facet_wrap(~organism, ncol=4, scales="free_y")+ 
  scale_fill_manual("Blast match type", values = rev(at_color_pals$hex_value[at_color_pals$colour == "N" & at_color_pals$intensity %in% c(800, 100)])) + 
  theme_figure + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(">=100 reads") + scale_x_discrete("Class") + scale_y_continuous("Count")


leg_5prime= ggplot_legend(com_inc_10)

pdf("plots/Supp_Figure_5primes.pdf", height=5,width=6.6)
ggdraw() +
  draw_plot(com_inc_10 + theme(legend.position = "none"), 0,0.5,0.8,0.5) + 
  draw_plot(com_inc_100+ theme(legend.position = "none"), 0,0,0.8,0.5) + 
  draw_plot(leg_5prime, 0.8,0.25,0.2,0.5) +
  draw_plot_label(c("A","B"), size=12, x=c(0), y=c(1,0.5))   
dev.off()


#UTRS


aggregate(read_counts ~ base_class, trinity_translated[trinity_translated$read_counts >=10,], summary)


trinity_translated$count_class = NA
trinity_translated$count_class[trinity_translated$read_counts >= 10] = "1e1"
trinity_translated$count_class[trinity_translated$read_counts >= 100] = "1e2"
trinity_translated$count_class[trinity_translated$read_counts >= 1000] = "1e3"
trinity_translated$count_class[trinity_translated$read_counts >= 10000] = "1e4"
trinity_translated$count_class[trinity_translated$read_counts >= 100000] = "1e5"

#### Suppl. figure

# this should be + strand only... negative strands are likely incorrect
pdf("plots/Supp_Figure_read_counts.pdf", width=6.7, height = 3)
ggplot(trinity_translated[trinity_translated$read_counts >=10,], 
       aes(x=count_class, y=lcs_percent_of_ens, fill=organism)) + 
  geom_boxplot(outlier.shape = NA, size=0.25) + facet_wrap(~organism, ncol=4) + 
  theme_figure + 
  theme(legend.position = "none") + scale_fill_brewer(palette = "Set2") + 
  scale_x_discrete("Read counts (greater than)") + 
  scale_y_continuous("Percent coverage of Ensembl ORF")
dev.off()

# Assembled 5' UTR length
ens_up = read.delim("data/ensembl_upstream_stops.txt")
ens_up$has_upstream_stop_f = ifelse(ens_up$has_upstream_stop == 1, T,F)

ens_up$organism = "h. sapiens"
ens_up$organism[ens_up$org == "zebrafish"] = "d. rerio"
ens_up$organism[ens_up$org == "arabidopsis"] = "a. thaliana"

# median 116 nt longer 5'utr
trinity_translated$has_upstream_stop_f = ifelse(trinity_translated$has_upstream_stop == 1, T,F)

get_accuracy = function(read_counts= 10){
  
  acc_df_all = NULL
  
  for(o in c("human", "arabidopsis", "yeast","zebrafish")){
    
    tt_get_cut = trinity_translated[trinity_translated$lcs_strand == "+" & 
                                      trinity_translated$org == o &
                                      trinity_translated$read_counts >= read_counts &
                                      trinity_translated$has_upstream_stop == 0 &
                                      trinity_translated$base_class_tag == "full_match" &
                                      trinity_translated$base_class %in% c("complete", "incomplete_3prime", "incomplete", "incomplete_5prime"),]
    
    tt_get_cut$base_class_5 = NA
    tt_get_cut$base_class_5[tt_get_cut$base_class %in% c("complete", "incomplete_3prime")] = "com"
    tt_get_cut$base_class_5[tt_get_cut$base_class %in% c("incomplete", "incomplete_5prime")] = "5ic"
    
    library(caret)
    acc = vector()
    accb = vector()
    for(i in 1:300){
      
      tt_get_cut$new_class = NA
      tt_get_cut$new_class[tt_get_cut$has_upstream_stop == 0 & tt_get_cut$base_class %in% c("complete", "incomplete_3prime") & tt_get_cut$upstream_stop > i] = "5ic"
      tt_get_cut$new_class[tt_get_cut$has_upstream_stop == 0 & tt_get_cut$base_class %in% c("complete", "incomplete_3prime") & tt_get_cut$upstream_stop <= i] = "com"
      tt_get_cut$new_class[tt_get_cut$base_class %in% c("incomplete", "incomplete_5prime") & tt_get_cut$first_M > i] = "5ic"
      tt_get_cut$new_class[tt_get_cut$base_class %in% c("incomplete", "incomplete_5prime") & tt_get_cut$first_M < i] = "com"
      
      cm  = confusionMatrix(factor(tt_get_cut$new_class), factor(tt_get_cut$base_class_5))
      acc[i] = cm$overall[1]
      accb[i] = cm$byClass[11]
    }
    
    acc_df = data.frame(org=o, acc,accb,upstream_len=1:300)
    acc_df_all = rbind(acc_df_all, acc_df)
  }
  
  
  acc_df_all$organism = "h. sapiens"
  acc_df_all$organism[acc_df_all$org == "zebrafish"] = "d. rerio"
  acc_df_all$organism[acc_df_all$org == "arabidopsis"] = "a. thaliana"
  acc_df_all$organism[acc_df_all$org == "yeast"] = "s. cerevisiae"
  
  return(acc_df_all)
  
}


acc_df_all_10 = get_accuracy(read_counts = 10)
acc_df_all_100 = get_accuracy(read_counts = 100)


get_max_cutoff = function(acc_df_all, metric = "acc"){
  ## also need ensembl upstream AAs // mean utr lengths?
  max_cutoff = data.frame(organism = c("h. sapiens","d. rerio","a. thaliana","s. cerevisiae"), upstream_len = NA)
  if(metric == "accb"){
    max_cutoff$upstream_len[1] = acc_df_all$upstream_len[acc_df_all$org=='human'][which.max(acc_df_all$accb[acc_df_all$org=='human'])]
    max_cutoff$upstream_len[2] = acc_df_all$upstream_len[acc_df_all$org=='zebrafish'][which.max(acc_df_all$accb[acc_df_all$org=='zebrafish'])]
    max_cutoff$upstream_len[3] = acc_df_all$upstream_len[acc_df_all$org=='arabidopsis'][which.max(acc_df_all$accb[acc_df_all$org=='arabidopsis'])]
    max_cutoff$upstream_len[4] = acc_df_all$upstream_len[acc_df_all$org=='yeast'][which.max(acc_df_all$accb[acc_df_all$org=='yeast'])]
    max_cutoff$accb = aggregate(accb~org, acc_df_all, max)[c(2,4,1,3),2]
  }else if(metric == "acc"){
    ## also need ensembl upstream AAs // mean utr lengths?
    max_cutoff = data.frame(organism = c("h. sapiens","d. rerio","a. thaliana","s. cerevisiae"), upstream_len = NA)
    max_cutoff$upstream_len[1] = acc_df_all$upstream_len[acc_df_all$org=='human'][which.max(acc_df_all$acc[acc_df_all$org=='human'])]
    max_cutoff$upstream_len[2] = acc_df_all$upstream_len[acc_df_all$org=='zebrafish'][which.max(acc_df_all$acc[acc_df_all$org=='zebrafish'])]
    max_cutoff$upstream_len[3] = acc_df_all$upstream_len[acc_df_all$org=='arabidopsis'][which.max(acc_df_all$acc[acc_df_all$org=='arabidopsis'])]
    max_cutoff$upstream_len[4] = acc_df_all$upstream_len[acc_df_all$org=='yeast'][which.max(acc_df_all$acc[acc_df_all$org=='yeast'])]
    max_cutoff$acc = aggregate(acc~org, acc_df_all, max)[c(2,4,1,3),2]
  }
  return(max_cutoff)
}

max_cutoff_10_acc = get_max_cutoff(acc_df_all_10, metric="acc")
max_cutoff_10_acb = get_max_cutoff(acc_df_all_10, metric="accb")
max_cutoff_100_acc = get_max_cutoff(acc_df_all_100, metric="acc")
max_cutoff_100_acb = get_max_cutoff(acc_df_all_100, metric="accb")


plot_accb_100 = 
  ggplot(acc_df_all_100, aes(x=upstream_len, y=accb, col=organism)) + 
  geom_line() + 
  theme_figure + 
  scale_x_continuous("Upstream length cutoff (AA)") +
  scale_y_continuous("Balanced accuracy") +
  scale_color_brewer(palette = "Set2", "Species") + 
  geom_segment(data=max_cutoff_100_acb, aes(x=upstream_len,y=0,xend=upstream_len,yend=accb), linetype=2) + 
  coord_cartesian(ylim = c(min(acc_df_all_100$accb), max(acc_df_all_100$accb)), xlim = c(0,100)) + 
  ggtitle("Balanced Accuracy\n(>=100 read counts)")

plot_acc_100 = 
  ggplot(acc_df_all_100, aes(x=upstream_len, y=acc, col=organism)) + 
  geom_line() + 
  theme_figure + 
  scale_x_continuous("Upstream length cutoff (AA)") +
  scale_y_continuous("Accuracy") +
  scale_color_brewer(palette = "Set2", "Species") + 
  geom_segment(data=max_cutoff_100_acc, aes(x=upstream_len,y=0,xend=upstream_len,yend=acc), linetype=2) + 
  theme(legend.position = "bottom") +
  guides(col=guide_legend(ncol=2)) +
  coord_cartesian(ylim = c(min(acc_df_all_100$acc), max(acc_df_all_100$acc)), xlim = c(0,200))+ ggtitle("Accuracy\n(complete vs. 5' incomplete)")

plot_accb_10 = 
  ggplot(acc_df_all_10, aes(x=upstream_len, y=accb, col=organism)) + 
  geom_line() + 
  theme_figure + 
  scale_x_continuous("Upstream length cutoff (AA)") +
  scale_y_continuous("Balanced accuracy") +
  scale_color_brewer(palette = "Set2", "Species") + 
  geom_segment(data=max_cutoff_10_acb, aes(x=upstream_len,y=0,xend=upstream_len,yend=accb), linetype=2) + 
  coord_cartesian(ylim = c(min(acc_df_all_10$accb), max(acc_df_all_10$accb)), xlim = c(0,100)) + 
  ggtitle("Balanced Accuracy\n(>=10 read counts)")

plot_acc_10 = 
  ggplot(acc_df_all_10, aes(x=upstream_len, y=acc, col=organism)) + 
  geom_line() + 
  theme_figure + 
  scale_x_continuous("Upstream length cutoff (AA)") +
  scale_y_continuous("Accuracy") +
  scale_color_brewer(palette = "Set2", "Species") + 
  geom_segment(data=max_cutoff_10_acc, aes(x=upstream_len,y=0,xend=upstream_len,yend=acc), linetype=2) + 
  guides(col=guide_legend(ncol=2)) +
  coord_cartesian(ylim = c(min(acc_df_all_10$acc), max(acc_df_all_10$acc)), xlim = c(0,200))+ ggtitle("Accuracy\n(>=10 read counts)") +
  guides(colour = guide_legend(nrow = 1)) + theme(legend.title.align=0.5)

### accuracy cutoff plots
leg_1 = ggplot_legend(plot_acc_10)

pdf("plots/Supp_Figure_cutoffs.pdf", height=3,width=6.6)
ggdraw() +
  draw_plot(plot_accb_100 + theme(legend.position = "none"), 0,0.1,0.33,0.9) + 
  draw_plot(plot_acc_10+ theme(legend.position = "none"), 0.33,0.1,0.33,0.9) + 
  draw_plot(plot_accb_10+ theme(legend.position = "none"), 0.66,0.1,0.33,0.9) + 
  draw_plot(leg_1, 0.25, 0, 0.5, 0.1) +
  draw_plot_label(c("A","B","C"), size=12, x=c(0,0.33, 0.66), y=c(1))   
dev.off()

#SUPP TABLE DATA
aggregate(upstream_stop_nt~  has_upstream_stop_f+organism, ens_up, median)%>% arrange(desc(organism), desc(has_upstream_stop_f))
aggregate(utr5_len_nt~ has_upstream_stop_f+organism, ens_up, median) %>% arrange(desc(organism), desc(has_upstream_stop_f))

dist_summary = aggregate((upstream_stop*3)+0.1 ~ organism+has_upstream_stop_f, ens_up, median)
colnames(dist_summary)[3] = "median"
dist_summary$y_max = 6800
dist_summary$y_max[dist_summary$has_upstream_stop_f == T] = 6500
dist_summary$hjust = ifelse(dist_summary$has_upstream_stop_f, 1,0)


plot_ensembl_upstream = 
  ggplot(ens_up, aes(fill=has_upstream_stop_f, x=(upstream_stop*3)+0.1)) + 
geom_histogram(position = "dodge", bins=30) + scale_x_log10("Distance to upstream STOP /\n5'UTR length (nt)") + 
  theme_figure + facet_wrap(~organism, ncol=2) +
  geom_vline(data=dist_summary, aes(xintercept=median, col=has_upstream_stop_f), linetype=2) +
  scale_fill_manual(values = c('#1f78b4','#b2df8a'), "STOP\nupstream\nof ORF") + 
  scale_color_manual(values = c('#1f78b4','#b2df8a'), "STOP\nupstream\nof ORF") +
  geom_label(data=dist_summary,aes(x=median, label=round(median, 0), y = y_max, fill=has_upstream_stop_f, hjust=hjust), alpha=0.2, size=2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.75, 0.25)) + coord_cartesian(xlim=c(1,1000), ylim=c(0,6800*1.1))


filter = which(trinity_translated$base_class_tag == "full_match" &
                 trinity_translated$base_class %in% c("complete", "incomplete_3prime") & 
                 trinity_translated$read_counts >= 100)

dist_summary_tt = aggregate((upstream_stop*3)+0.1 ~ organism+has_upstream_stop_f, trinity_translated[filter,], median)
colnames(dist_summary_tt)[3] = "median"
dist_summary_tt$y_max = 2700
dist_summary_tt$y_max[dist_summary_tt$has_upstream_stop_f == T] = 2600
dist_summary_tt$hjust = ifelse(dist_summary_tt$has_upstream_stop_f, 1,0)

plot_ass_up = ggplot(trinity_translated[filter,],
       aes(fill=has_upstream_stop_f, x=(upstream_stop*3)+0.1)) + 
  geom_histogram(position = "dodge", bins=30) + scale_x_log10("Distance to upstream STOP /\n5'UTR length (nt)") + 
  theme_figure + facet_wrap(~organism, ncol=2) +
  geom_vline(data=dist_summary_tt, aes(xintercept=median, col=has_upstream_stop_f), linetype=2) +
  scale_fill_manual(values = c('#1f78b4','#b2df8a'), "STOP upstream of ORF") + 
  scale_color_manual(values = c('#1f78b4','#b2df8a'), "STOP upstream of ORF") +
  geom_label(data=dist_summary_tt,aes(x=median, label=round(median, 0), y = y_max, fill=has_upstream_stop_f, hjust=hjust), alpha=0.2, size=2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_cartesian(xlim=c(1,1000), ylim=c(0,2700*1.1))



aggregate(upstream_stop*3 ~ has_upstream_stop_f+organism, trinity_translated[filter,], median) %>% 
  arrange(desc(organism), desc(has_upstream_stop_f))
aggregate(lcs_start_site_nt ~ has_upstream_stop_f+organism, trinity_translated[filter,], median) %>% 
  arrange(desc(organism), desc(has_upstream_stop_f))

table(trinity_translated$organism[filter], trinity_translated$has_upstream_stop_f[filter]) %>% 
  as.data.frame() %>% arrange(desc(Var1), desc(Var2))
table(trinity_translated$organism[filter], trinity_translated$has_upstream_stop_f[filter]) / rowSums(table(trinity_translated$organism[filter], trinity_translated$has_upstream_stop_f[filter]))

aggregate(lcs_start_site_nt ~ has_upstream_stop_f+organism, trinity_translated[filter,], median)

plot_5utr_len = ggplot(trinity_translated[filter,], 
                       aes(y=lcs_start_site_nt, x=organism, fill=organism)) + 
  geom_boxplot(outlier.shape = NA, lwd=0.25) + 
  scale_y_log10("Assembled 5'UTR length (nt)") + 
  scale_fill_brewer(palette = "Set2") + 
  scale_x_discrete("Species") +
  theme_figure + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  easy_remove_legend() + ggtitle("Assembled 5'UTR lengths")

plot_5utr_diffs= ggplot(trinity_translated[filter,], 
                        aes(y=lcs_start_site_nt-ens_5utr_len, x=organism,fill=organism)) + 
  geom_boxplot(outlier.shape = NA, lwd=0.25)+ 
  coord_cartesian(ylim = c(-1500, 1500)) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_figure + easy_remove_legend()+ 
  scale_x_discrete("Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous("Assembled 5'UTR -\nreference 5'UTR (nt)")+ ggtitle("5' UTR length\ndifferences")

pdf("plots/Figure3_additional_plots.pdf", height=5, width=4.4)
ggdraw() +
  draw_plot(plot_ensembl_upstream + ggtitle("Reference transcripts"), 0,0.4,0.5,0.6) + 
  draw_plot(plot_ass_up + theme(legend.position = "none")+ ggtitle("Assembled transcripts"), 0.5,0.4,0.5,0.6) + 
  draw_plot(plot_acc_100+ theme(legend.position = "bottom"), 0.4,0,0.6,0.4) +
  draw_plot(plot_5utr_diffs+ theme(legend.position = "none"), 0,0,0.4,0.4) +
  draw_plot_label(c("A","B","C","D"), size=12, x=c(0,0.5,0,0.4), y=c(1,1,0.4,0.4)) 
dev.off()  

pdf("plots/Figure3.pdf", height=2.5, width=3.3)
plot_acc_100 + theme(legend.position = "right")+
  guides(col=guide_legend(ncol=1))
dev.off()  


 ############################

# BORF vs. Transdecoder vs. known

#write files for 'expressed' protein seqs
# for(org in c('yeast',"human", "arabidopsis", 'zebrafish')){
# 
#     tt = fread(paste0('data/', org, "/trinity_translated_both_strands.txt"), data.table = F)
# 
#     borf_results = fread(paste0("data/", org, "/borf/Trinity.txt"), data.table = F)
#     borf_seq = halpme::read_fasta2df(paste0("data/", org, "/borf/Trinity.pep"))
# 
#     borf_results = cbind(borf_results, borf_seq[,2])
#     borf_results = borf_results[borf_results$transcript_id %in% tt$subject,]
#     borf_results = cbind(borf_results, tt[match(borf_results$transcript_id, tt$subject),])
# 
#     borf_results = borf_results[borf_results$read_counts >=100,]
# 
#     borf_seq_out = borf_results[,c(1,15)]
#     write.table(borf_seq_out, file=paste0("data/", org, "/borf/borf_expressed.pep"), quote=F,row.names = F,col.names = F,sep="\n")
# 
# 
#     transdecoder_results = halpme::read_fasta2df(paste0("data/", org, "/Trinity.fasta.transdecoder_dir/longest_orfs.pep"))
#     transdecoder_results$orf_id = strv_split(transdecoder_results$seq_id, "[ ]", 1)
#     transdecoder_results$transcript_id = strv_split(transdecoder_results$orf_id, "[.]p", 1)
#     transdecoder_results$strand = strv_split2(transdecoder_results$seq_id, "[(]", "[)]")
# 
#     locs = strv_split(transdecoder_results$seq_id, "[ ]", 5)
#     transdecoder_results$start_site_nt = as.numeric(strv_split2(locs, "[:]", "[-]"))
#     transdecoder_results$stop_site_nt = as.numeric(strv_split2(locs, "[-]", "[(]"))
#     transdecoder_results$orf_length_aa = as.numeric(strv_split2(transdecoder_results$seq_id, "len[:]", "[ ]"))
#     transdecoder_results$orf_length_aa = nchar(gsub("[*]", "", transdecoder_results$seq))
# 
#     transdecoder_results$orf_class = strv_split2(transdecoder_results$seq_id, "type[:]", "[ ]")
# 
#     transdecoder_results = transdecoder_results[transdecoder_results$transcript_id %in% tt$subject,]
#     transdecoder_results = cbind(transdecoder_results, tt[match(transdecoder_results$transcript_id, tt$subject),])
#     transdecoder_results = transdecoder_results[transdecoder_results$read_counts >=100,]
# 
#     transdecoder_seqs=transdecoder_results[,c(1,2)]
#     transdecoder_seqs$seq_id = paste0(">", transdecoder_seqs$seq_id)
#     transdecoder_seqs$seq = gsub("[*]","",transdecoder_seqs$seq)
#     write.table(transdecoder_seqs, file=paste0("data/", org, "/Trinity.fasta.transdecoder_dir/transdecoder_expressed.pep"), quote=F,row.names = F,col.names = F,sep="\n")
# }

org = "human"
tt_all=NULL
td_all=NULL
for(org in c("human", "arabidopsis", 'zebrafish', 'yeast')){
  
  tt = fread(paste0('data/', org, "/trinity_translated_both_strands.txt"), data.table = F)
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
  write.table(seqs_b_pos, file=paste0("data/",org, "/borf_interpro_pos.fa"), quote=F, row.names = F, col.names = F, sep='\n')
  
  
  
  
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
  
  write.table(seqs_td_pos, file=paste0("data/",org, "/transdecoder_interpro_pos.fa"), quote=F, row.names = F, col.names = F, sep='\n')
  write.table(seqs_td_neg, file=paste0("data/",org, "/transdecoder_interpro_neg.fa"), quote=F, row.names = F, col.names = F, sep='\n')
  
  
  #transdecoder_results = transdecoder_results[transdecoder_results$orf_length_aa >=100,]
  
  
  colnames(borf_results) = paste0(colnames(borf_results), "_borf")
  #tt = left_join(tt, borf_results, by=c('subject'='transcript_id_borf'))
  tt = cbind(tt, borf_results[match(tt$subject, borf_results$transcript_id_borf),])
  
  
  colnames(transdecoder_results) = paste0(colnames(transdecoder_results), "_transdecoder")
  #tt = left_join(tt, transdecoder_results[,], by=c('subject'='transcript_id_transdecoder'))
  tt = cbind(tt, transdecoder_results[match(tt$subject, transdecoder_results$transcript_id_transdecoder),])
  
  tt$orf_class_transdecoder[is.na(tt$orf_class_transdecoder)] = "no_ORF"
  tt$orf_class_transdecoder[which(tt$orf_class_transdecoder == "3prime_partial")] = "incomplete_3prime"
  tt$orf_class_transdecoder[which(tt$orf_class_transdecoder == "5prime_partial")] = "incomplete_5prime"
  tt$orf_class_transdecoder[which(tt$orf_class_transdecoder == "internal")] = "incomplete"
  
  tt$orf_class_borf[is.na(tt$orf_class_borf)] = "no_ORF"
  tt$org=org
  tt_all=plyr::rbind.fill(tt_all, tt)
  
}

tt_all$organism = "h. sapiens"
tt_all$organism[tt_all$org == "zebrafish"] = "d. rerio"
tt_all$organism[tt_all$org == "arabidopsis"] = "a. thaliana"
tt_all$organism[tt_all$org == "yeast"] = "s. cerevisiae"

td_all$organism = "h. sapiens"
td_all$organism[td_all$org == "zebrafish"] = "d. rerio"
td_all$organism[td_all$org == "arabidopsis"] = "a. thaliana"
td_all$organism[td_all$org == "yeast"] = "s. cerevisiae"


tt_all$strand_transdecoder[which(is.na(tt_all$strand_transdecoder))] = "+"
tt_all$strand_transdecoder[which((tt_all$strand_transdecoder == "*"))] = "+"

td_all$strand[which(is.na(td_all$strand))] = "+"
td_all$strand[which((td_all$strand == "*"))] = "+"

tt_all_exp = tt_all[tt_all$read_counts >= 100,]
td_all_exp = td_all[td_all$read_counts >= 100,]


# table of transdecoder 99AA ORFs
table(td_all_exp$org, td_all_exp$orf_length_aa < 100)
table(td_all_exp$org, td_all_exp$orf_length_aa < 100) / rowSums(table(td_all_exp$org, td_all_exp$orf_length_aa < 100))

####
td_all = arrange(td_all, org, desc(orf_length_aa))
td_all_exp = td_all[td_all$read_counts>= 100,]
td_all_exp$longest_orf = ifelse(duplicated(paste0(td_all_exp$org, td_all_exp$transcript_id)), "Shorter ORF", "Longest ORF")
td_all_exp$longer_short_blastpmatch = paste0(td_all_exp$longest_orf, " - ", ifelse(td_all_exp$base_in_blastp == 1 | td_all_exp$base_in_blastp_gene == 1, "Correct BlastP match", "Incorrect / no match"))
td_all_exp$base_in_blastp_gene_no = ifelse(td_all_exp$base_in_blastp_gene == 1, "correct hit", "incorrect hit")
td_all_exp$base_in_blastp_gene_no[which(is.na(td_all_exp$hit_blastp))] = "no hit"


transdecoder_table = as.data.frame(table(td_all_exp$organism))
colnames(transdecoder_table)= c("Species", "total_orfs")
transdecoder_table$total_transcripts = as.data.frame(table(td_all_exp$organism[!duplicated(paste0(td_all_exp$organism, td_all_exp$transcript_id))]))[,2]

# multiple ORFs per tx
transdecoder_table = td_all_exp %>% filter(longest_orf == "Shorter ORF") %>% filter(!duplicated(paste0(org, transcript_id))) %>% select(organism) %>% table() %>% as.data.frame() %>% 
  rename(c('.'='Species', 'Freq'='total_transcripts_multi_orf')) %>% left_join(transdecoder_table,.,by="Species")
# txs with longest orf on - strand
transdecoder_table = td_all_exp %>% filter(longest_orf == "Longest ORF" & strand=="-") %>% select(organism) %>% table() %>% as.data.frame() %>% 
  rename(c('.'='Species', 'Freq'='longest_ORF_neg')) %>% left_join(transdecoder_table,.,by="Species")
# txs with incorrect longest orf on - strand
transdecoder_table = td_all_exp %>% filter(longest_orf == "Longest ORF" & strand=="-"& base_in_blastp_gene_no != "correct hit") %>% select(organism) %>% table() %>% as.data.frame() %>% 
  rename(c('.'='Species', 'Freq'='longest_ORF_neg_incorrect')) %>% left_join(transdecoder_table,.,by="Species")

transdecoder_table$incorrect_neg_percent = (transdecoder_table$longest_ORF_neg_incorrect / transdecoder_table$total_transcripts)*100


# multiple ORFs per tx
transdecoder_table = td_all_exp %>% filter(base_in_blastp_gene_no == "incorrect hit") %>% select(organism) %>% table() %>% as.data.frame() %>% 
  rename(c('.'='Species', 'Freq'='false_pos_hits')) %>% left_join(transdecoder_table,.,by="Species")
transdecoder_table = td_all_exp %>% filter(base_in_blastp_gene_no == "no hit") %>% select(organism) %>% table() %>% as.data.frame() %>% 
  rename(c('.'='Species', 'Freq'='no_hits')) %>% left_join(transdecoder_table,.,by="Species")

transdecoder_table$increase_blast_search = transdecoder_table$total_orfs / transdecoder_table$total_transcripts
transdecoder_table$false_pos_percent = (transdecoder_table$false_pos_hits / transdecoder_table$total_orfs)*100
transdecoder_table$no_hit_percent = (transdecoder_table$no_hits / transdecoder_table$total_orfs)*100

transdecoder_table.longest = transdecoder_table[,c(1,3,5,6,7)]
write_csv(transdecoder_table.longest, "data/table_s_transdecoder_longest.csv")
transdecoder_table.all = transdecoder_table[,c(1,3,2,4,8,9,11,12,10)]
write_csv(transdecoder_table.all, "data/table_s_transdecoder_all.csv")


table(td_all_exp$organism)
table(td_all_exp$organism[!duplicated(paste0(td_all_exp$transcript_id, td_all_exp$organism))])

td_all_shorter = td_all_exp[td_all_exp$longest_orf == "Shorter ORF",]
td_all_shorter = td_all_shorter[!duplicated(paste0(td_all_shorter$transcript_id, td_all_shorter$organism)),]
table(td_all_shorter$organism)

td_strand_plot_data = as.data.frame(table(td_all_exp$longer_short_blastpmatch, td_all_exp$strand, td_all_exp$organism))

td_strand_plot_data$Freq[td_strand_plot_data$Var2 == "-"] = (td_strand_plot_data$Freq *-1)[td_strand_plot_data$Var2 == "-"]                                       

td_strand_plot_data %>% filter(Var1 == "Longest ORF - Incorrect / no match") %>% 
  filter(Var2 == "-")

text_df = data.frame(text = c("Sense (+) ORFs", "Antisense (-) ORFs"), 
                     Var1 = rep('Longest ORF - Incorrect / no match',2),
                     Var2 = c("+", "-"),
                     y = c(20000, -20000),
                     Var3 = rep("a. thaliana",2))

plot_strand_orfs = 
  ggplot(td_strand_plot_data, aes(x= Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity", position = "identity") + facet_wrap(~Var3, scales = "free_x", ncol = 1) + geom_hline(yintercept = 0) + 
  geom_blank(aes(y = -Freq)) + 
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_manual(values = c('#6a3d9a','#e31a1c', '#cab2d6','#fb9a99'),"") + 
  theme_figure +
  scale_x_discrete("Transdecoder ORF", 
                   labels = c("Longest/Correct","Longest/Incorrect",
                              "Shortest/Correct","Shortest/Incorrect")) + 
  scale_y_continuous("Number of ORFs") + 
  geom_text(data = text_df, aes(x=Var1, y=y, label = text), size=2) + 
  coord_flip() + ggtitle("Transdecoder all ORF BlastP matches")+ 
  theme(panel.spacing = unit(0.1, "lines"))

pdf("plots/Figure4.pdf", height=3.5,width=5.5)
plot_strand_orfs
dev.off()


aggregate(Freq ~ Var3, td_strand_plot_data, function(x) sum(abs(x)))

# longest only - incorrect gene match
longest_match = as.data.frame(table(td_all_exp$base_in_blastp_gene[td_all_exp$longest_orf == "Longest ORF"], 
                                    td_all_exp$strand[td_all_exp$longest_orf == "Longest ORF"], 
                                    td_all_exp$organism[td_all_exp$longest_orf == "Longest ORF"]))
# incorrect/no gene match
td_all_exp$base_in_blastp_gene_no = ifelse(td_all_exp$base_in_blastp_gene == 1, "correct hit", "incorrect hit")
td_all_exp$base_in_blastp_gene_no[which(is.na(td_all_exp$hit_blastp))] = "no hit"

table(td_all_exp$organism[which(td_all_exp$longest_orf == "Longest ORF" & 
                                  td_all_exp$strand =="-")])

as.data.frame(table(td_all_exp$base_in_blastp_gene_no== 'correct hit', 
                     td_all_exp$organism))

all_match = as.data.frame(table(td_all_exp$base_in_blastp_gene_no, 
                                td_all_exp$strand, 
                                td_all_exp$organism))

corr_incor_no = as.data.frame(table(td_all_exp$base_in_blastp_gene_no, 
                                    td_all_exp$strand, td_all_exp$organism))


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
  scale_x_discrete("Species") + scale_y_continuous("Accuracy")

colnames(total_ts)[1:2] = c("Species/Annotation", "Total transcripts")
write_csv(total_ts, "supp_compare_accuracy.csv")

runtimes = read.delim("data/timing_tests/timing_tests.txt")
plot_runtime = 
  ggplot(runtimes, aes(x=n, y=time_sec, col=software)) + 
  geom_point() +  scale_color_manual(values = c("#00B8D9", "#FFAB00"), "Method")+ 
  #scale_x_log10() + scale_y_log10() + 
  geom_smooth() +
  #geom_point(data = runtimes[runtimes$set != "random",], aes(shape=set), alpha=1, size=2) + 
  theme_figure + 
  scale_x_continuous("Transcripts") + scale_y_continuous("Total time (seconds)")+ 
  ggtitle("Run times")


leg_tb = ggplot_legend(plot_method_accuracy)

pdf("plots/Figure5.pdf", height=5, width=4.4)
ggdraw() +
  draw_plot(plot_corrincorr + theme(panel.spacing.y=unit(0, "lines"), legend.position = "right") + guides(fill=guide_legend(ncol=1)), 0,0.4,1,0.6) + 
  draw_plot(plot_method_accuracy + theme(legend.position ="none", axis.text.x = element_text(angle = 45, hjust = 1)), 0,0,0.4,0.4) +
  draw_plot(plot_runtime + theme(legend.position ="none",axis.text.x = element_text(angle = 45, hjust = 1)), 0.4,0,0.4,0.4) +
  draw_plot(leg_tb, 0.8,0,0.2,0.4) +
  draw_plot_label(c("A","B", "C"), size=12, x=c(0,0,0.4), y=c(1,0.4,0.4)) 
dev.off()



## tt_all
write_delim(tt_all, "data/trinity_translated_all.txt", delim = "\t")