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
library(ggpubr)

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

tt_yeast = fread("data/yeast/trinity_translated_both_strands.txt", data.table = F)
tt_ara = fread("data/arabidopsis/trinity_translated_both_strands.txt", data.table = F)
tt_zebrafish = fread("data/zebrafish/trinity_translated_both_strands.txt", data.table = F)
tt_human = fread("data/human/trinity_translated_both_strands.txt", data.table = F)

tt_yeast$org = "yeast"
tt_ara$org = "arabidopsis"
tt_zebrafish$org = "zebrafish"
tt_human$org = "human"

coding_potential = fread("data/yeast/rnasamba/rnasamba_partial.tsv", data.table = F)
tt_yeast$coding_potential = coding_potential$coding_score[match(tt_yeast$subject, halpme::strv_split(coding_potential$sequence_name, "[ ]", 1))]
coding_potential = fread("data/arabidopsis/rnasamba/rnasamba_partial.tsv", data.table = F)
tt_ara$coding_potential = coding_potential$coding_score[match(tt_ara$subject, halpme::strv_split(coding_potential$sequence_name, "[ ]", 1))]
coding_potential = fread("data/zebrafish/rnasamba/rnasamba_partial.tsv", data.table = F)
tt_zebrafish$coding_potential = coding_potential$coding_score[match(tt_zebrafish$subject, halpme::strv_split(coding_potential$sequence_name, "[ ]", 1))]
coding_potential = fread("data/human/rnasamba/rnasamba_partial.tsv", data.table = F)
tt_human$coding_potential = coding_potential$coding_score[match(tt_human$subject, halpme::strv_split(coding_potential$sequence_name, "[ ]", 1))]



trinity_translated_bs = plyr::rbind.fill(tt_yeast[tt_yeast$read_counts >=100,],tt_ara,tt_zebrafish,tt_human)

trinity_translated_bs$base_class_factor = factor(trinity_translated_bs$base_class, levels = rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage")))


trinity_translated_bs$organism = "h. sapiens"
trinity_translated_bs$organism[trinity_translated_bs$org == "zebrafish"] = "d. rerio"
trinity_translated_bs$organism[trinity_translated_bs$org == "arabidopsis"] = "a. thaliana"
trinity_translated_bs$organism[trinity_translated_bs$org == "yeast"] = "s. cerevisiae"

at_color_pals = read.delim("at_colour_pal.txt")
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

trinity_translated_bs$base_class_f2 = factor(gsub("_", " ", trinity_translated_bs$base_class), levels = c(levels = gsub("_"," ",rev(c('complete', "incomplete_5prime","incomplete_3prime", "incomplete", 'complete_partial', "incomplete_5prime_partial","incomplete_3prime_partial", "incomplete_partial", "insufficient_coverage")))))


trinity_translated_bs$match_strand = factor(trinity_translated_bs$match_strand, levels = c("+",  "-"))
s1 = 
    ggplot(trinity_translated_bs[which(trinity_translated_bs$read_counts >=100),], aes(x=match_strand, fill=base_class_f2)) + 
    geom_bar() + facet_wrap(~organism, scales="free", ncol=4) +
    theme_figure +
    scale_x_discrete("BlastN match strand") + 
    scale_y_continuous("Number of transcripts") +
    theme(legend.position = "bottom") +
    scale_fill_manual(values=base_class_pal2, "Assembled transcript class") + ggtitle("ORF match type") + 
        guides(fill=guide_legend(title.position = "top", title.hjust = 0.5, nrow=2))

ggplot(trinity_translated_bs, aes(fill=match_strand, x=read_counts+0.01)) + 
    geom_histogram(bins=50, position='dodge') + 
    facet_wrap(~org, scales="free", ncol=4) + 
    theme_figure + scale_x_log10()



s2 = ggplot(trinity_translated_bs[which(trinity_translated_bs$read_counts >=100),], 
       aes(x=match_strand, y=read_counts+0.01)) + 
    geom_boxplot( lwd=0.25, aes(fill=match_strand), outlier.shape = NA) + facet_wrap(~organism, scales="free", ncol=4) + 
    scale_x_discrete("BlastN match strand") + scale_fill_manual(values = c("#1f78b4","#a6cee3")) +
    scale_y_log10("Total read counts",expand = expansion(mult = c(0.05, 0.15))) +
    stat_compare_means(comparisons = list(c("+", "-")), method = "wilcox.test", size = 2, tip.length = 0.01, vjust=-0.2, bracket.size = 0.25) +
    theme_figure + theme(legend.position = "none") + ggtitle("Expression")






s3 = ggplot(trinity_translated_bs[which(trinity_translated_bs$base_class != "insufficient_coverage" & trinity_translated_bs$read_counts >=100),], 
            aes(fill=match_strand, y=coding_potential, x=match_strand)) + 
    scale_fill_manual(values =  c("#1f78b4","#a6cee3")) +
    scale_y_continuous("RNAsamba coding score",expand = expansion(mult = c(0.05, 0.15))) + 
    scale_x_discrete("BlastN match strand") +
    geom_boxplot(outlier.shape = NA, lwd=0.25) + facet_wrap(~organism, scales="free", ncol=4) + 
    theme_figure + theme(legend.position = "none") + ggtitle("Coding probability") + 
    stat_compare_means(comparisons = list(c("+", "-")), method = "wilcox.test", size = 2, tip.length = 0.01, vjust=-0.2, bracket.size = 0.25)


library(cowplot)

pdf("plots/Figure2.pdf", height=5.5,width=4.4)
ggdraw() +
    draw_plot(s1, 0,0.6,1,0.4) + 
    draw_plot(s2 + theme(legend.position = "none"), 0,0.3,1,0.3)+
    draw_plot(s3 + theme(legend.position = "none"), 0,0,1,0.3)+
    #draw_plot(leg, 0.86,0,0.14,1) +
    draw_plot_label(c("A","B", "C"), size=12, x=c(0,0.,0), y=c(1,0.6,0.3))   
dev.off()
