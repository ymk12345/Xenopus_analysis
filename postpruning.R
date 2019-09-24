.libPaths("/moto/ziab/users/yk2840/software/R_Libs/")

setwd("/moto/ziab/users/yk2840/Frog_workspace/HMM/Ancestry_nofilter_results/posterior/")

library(parallel)
numCores<-detectCores()

files<-list.files(pattern = "posterior$")

index <- read.table("/moto/ziab/users/yk2840/Frog_workspace/HMM/Ancestry_nofilter_results/posterior/snps_thinned_frac0.2.txt", 
                           quote="\"", comment.char="", stringsAsFactors=FALSE)[,1]
index<-index[order(index)]

library(foreach)

read_samples_pruned<-function(file1, index){
  df <- read.delim(file1, stringsAsFactors=FALSE)[index,]
  colnames(df)[3:5]<-c("L", "LP", "P")
  message("loaded")
  r<-mclapply(1:nrow(df),
              FUN=function(z){
                
                if(df[z,3]==df[z,4] & df[z,4]==df[z,5]){
                  gt<-"-"
                }else{
                  x<-which.max(df[z,3:5])
                  if(df[z,3:5][x]<.4){
                    gt<-"-"
                  }else{
                    gt<-colnames(df)[x+2]
                  }
                }
              },
              mc.cores = numCores)
  r<-unlist(r)
  r<-c(file1[1], r)
  return(r)
  
  
}


data2 <- foreach(i = files, .combine = rbind) %dopar% read_samples_pruned(i, index)

write.table(data2, paste0("data_thinned_filter_frac.2_072919.txt"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)


