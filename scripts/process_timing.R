library(tidyverse)
library(data.table)
library(halpme)

timing_out = fread(cmd = "grep real data/timing_tests/time.log", header = F, data.table = F)
timing_err = fread("data/timing_tests/command.out", header = F, sep='\n', data.table = F)

file_processed = timing_err[grepl("^Results in", timing_err$V1) | grepl(".fa$", timing_err$V1),]
file_processed = data.frame(file = file_processed, software=NA)
file_processed$file[grep("Results", file_processed$file)]  = strv_split2(file_processed$file[grep("Results", file_processed$file)], "and ", ".txt")

file_processed = file_processed[rep(1:nrow(file_processed), each=4),]
file_processed$software = rep(c("Borf", "Transdecoder","TransDecoder", "GeneMarkS-T"))
file_processed$sample_n = rep(seq(1, nrow(file_processed)/4), each=2)

file_processed$time = timing_out$V2
file_processed$time_sec = as.numeric(strv_split(file_processed$time, "m", 1)) * 60 + as.numeric(strv_split2(file_processed$time, "m", "s"))

# combine the two times from transdecoder
transdecoder_bothtimes = aggregate(time_sec~file, file_processed[file_processed$software == "Transdecoder",], sum)
file_processed = file_processed[!duplicated(paste0(file_processed$file, file_processed$software)),]
m = match(paste0(transdecoder_bothtimes$file, "Transdecoder"), paste0(file_processed$file, file_processed$software))
file_processed$time_sec[m] = transdecoder_bothtimes$time_sec

file_processed$n = as.numeric(strv_split2(file_processed$file, "_n", "_"))
file_processed = arrange(file_processed, sample_n)

write_delim(file_processed, "data/timing_tests/timing_tests.txt", delim='\t')
