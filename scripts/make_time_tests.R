library(halpme)
options(scipen =999)

human_trinity = read_fasta2df("data/human/Trinity.fasta")
n_txs = c(2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)

set.seed(1)

for(j in 1:5){
    for(i in seq_along(n_txs)){
        trinity_sample = human_trinity[sample(1:nrow(human_trinity), n_txs[i]),]
        trinity_sample$seq_id = paste0(">",trinity_sample$seq_id)
        write.table(trinity_sample, paste0("data/timing_tests/trinity_sample_n", n_txs[i], "_rep", j, ".fa"), quote=F, row.names = F,col.names = F, sep='\n')
    }
}

timing_script = vector()
for(i in seq_along(n_txs)){
    for(j in 1:5){
        
        timing_script = c(timing_script, "borf -h")
        timing_script = c(timing_script,"rm -rf *transdecoder_dir*")
        
        timing_script = c(timing_script,paste0("{ ", paste0("time borf -f ", "trinity_sample_n", n_txs[i], "_rep", j, ".fa"), " >> command.out ;} 2>> time.log"))
        
        timing_script = c(timing_script, paste0("echo trinity_sample_n", n_txs[i], "_rep", j, ".fa", " >> command.out"))
        
        timing_script = c(timing_script,paste0("{ ", paste0("time TransDecoder.LongOrfs -t ", "trinity_sample_n", n_txs[i], "_rep", j, ".fa"), " >> command.out ;} 2>> time.log"))

    }
}
    
write.table(timing_script, "data/timing_tests/time_orfs.sh", quote=F, row.names = F, col.names = F, sep='\n')
