library(ggplot2)

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

options(stringsAsFactors = F)


scholar_searches_denovo <- read.csv("scholar_searches_denovo.csv", header=TRUE)

scholar_searches_denovo[7,-1] / apply(scholar_searches_denovo[-7,-1], 2 ,  function(x) max(x, na.rm = T))

scholar_searches_denovo = 
scholar_searches_denovo %>% 
    pivot_longer(-X, names_to ="year", values_to="n_pubs") %>% 
    mutate(year = gsub("X", "", year))

pdf("plots/Supp_Figure_citations.pdf", width=4.5, height=2.5, useDingbats = F)
ggplot(scholar_searches_denovo, aes(x=year, col=X, y=n_pubs, group=X)) + 
    geom_line() + geom_point() + 
    scale_y_continuous("Number of citations") + scale_x_discrete("Year") + scale_color_discrete("Software Publication") + theme_figure
dev.off()

