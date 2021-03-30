library(tidyverse)

org_version = "human"
kallisto_files = list.files(paste0("data/", org_version, "/kallisto/"), pattern = ".kallisto")
samples = gsub("[.]kallisto", "", kallisto_files)

for(f in seq_along(kallisto_files)){
    
    kallisto_e = fread(paste0("data/", org_version, "/kallisto/", kallisto_files[f], "/abundance.tsv"))
    
    if(f==1){
        
        kallisto_est = kallisto_e[,c(1,4)]
    }else{
        kallisto_est = full_join(kallisto_est, kallisto_e[,c(1,4)], by='target_id')
    }
    
}
colnames(kallisto_est)[-1] = samples
write_tsv(kallisto_est, paste0("data/", org_version, "/kallisto/kallisto_est_counts.tsv"))
