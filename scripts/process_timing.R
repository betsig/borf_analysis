library(tidyverse)
library(data.table)
library(halpme)

timing_out = fread(cmd = "grep real data/timing_tests/time.log", header = F, data.table = F)
timing_err = fread("data/timing_tests/command.out", header = F, sep='\n', data.table = F)

file_processed = timing_err[grepl("^Results in", timing_err$V1) | grepl(".fa$", timing_err$V1),]
file_processed = data.frame(file = file_processed, software=NA)
file_processed$software[grep("Results", file_processed$file)] = "Borf"
file_processed$software[!grepl("Results", file_processed$file)] = "Transdecoder"
file_processed$file[grep("Results", file_processed$file)]  = strv_split2(file_processed$file[grep("Results", file_processed$file)], "and ", ".txt")
file_processed$file = gsub(".fa", "", file_processed$file)

file_processed$time = timing_out$V2
file_processed$time_sec = as.numeric(strv_split(file_processed$time, "m", 1)) * 60 + as.numeric(strv_split2(file_processed$time, "m", "s"))
file_processed$n = as.numeric(strv_split2(file_processed$file, "_n", "_"))

write_delim(file_processed, "data/timing_tests/timing_tests.txt", delim='\t')
