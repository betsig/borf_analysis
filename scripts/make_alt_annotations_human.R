library(GenomicRanges)
library(halpme)

gtf = rtracklayer::import("data/human/ensembl/Homo_sapiens.GRCh38.98.gtf")

## protein_coding only
pc_tx_ids = unique(paste0(gtf$transcript_id[which(gtf$type=="transcript" & gtf$transcript_biotype == "protein_coding")], ".", gtf$transcript_version[which(gtf$type=="transcript" & gtf$transcript_biotype == "protein_coding")]))
pc_pep_ids = unique(paste0(gtf$protein_id[which(gtf$type=="CDS" & gtf$transcript_biotype == "protein_coding")], ".", gtf$protein_version[which(gtf$type=="CDS" & gtf$transcript_biotype == "protein_coding")]))

# extract protein_coding transcripts from cdna.fa
write.table(pc_tx_ids, file="data/human/pc_tx_only_ids.txt", sep='\n', quote=F, row.names = F, col.names = F)
# extract protein_coding transcripts from .pep
write.table(pc_pep_ids, file="data/human/pc_pep_only_ids.txt", sep='\n', quote=F, row.names = F, col.names = F)

pc_tx_ids = unique(gtf$transcript_id[which(gtf$type=="transcript" & gtf$transcript_biotype == "protein_coding")])
pc_pep_ids = unique(gtf$protein_id[which(gtf$type=="transcript" & gtf$transcript_biotype == "protein_coding")])


pc_gene_ids = unique(gtf$gene_id[which(gtf$type=="transcript" & gtf$transcript_biotype == "protein_coding")])

gtf_pc = gtf[which((gtf$transcript_id %in% pc_tx_ids) | (gtf$type=="gene" & gtf$gene_id %in% pc_gene_ids))]
rtracklayer::export(gtf_pc, "data/human/human_pc.gtf", format = "GTF")


## CDS only
CDS_regions = as.data.frame(gtf_pc[which(gtf_pc$type == "CDS")])
CDS_tx_ranges = aggregate(start ~ transcript_id, CDS_regions, min)
CDS_tx_ranges$end = aggregate(end ~ transcript_id, CDS_regions, max)[,2]

CDS_gene_ranges = aggregate(start ~ gene_id, CDS_regions, min)
CDS_gene_ranges$end = aggregate(end ~ gene_id, CDS_regions, max)[,2]

gtf_cds = gtf_pc
gtf_cds = gtf_cds[(gtf_cds$type %in% c("gene", "transcript", "CDS"))]
gtf_cds$type[gtf_cds$type == "CDS"] = "exon"

mt = match(gtf_cds$transcript_id[gtf_cds$type == "transcript"], CDS_tx_ranges$transcript_id)
ranges(gtf_cds[gtf_cds$type == "transcript"]) = IRanges(start=CDS_tx_ranges$start[mt], end=CDS_tx_ranges$end[mt])

mg = match(gtf_cds$gene_id[gtf_cds$type == "gene"], CDS_gene_ranges$gene_id)
ranges(gtf_cds[gtf_cds$type == "gene"]) = IRanges(start=CDS_gene_ranges$start[mg], end=CDS_gene_ranges$end[mg])

rtracklayer::export(gtf_cds, "data/human/human_cds.gtf", format = "GTF")

# gtf2bed (from how_are_we_stranded here)
#gtf2bed --gtf human_cds.gtf --bed human_cds.bed
#bedtools getfasta -fi Homo_sapiens.GRCh38.dna.toplevel.fa -bed human_cds.bed -name -split -s -fo human_cds.cdna.fa



