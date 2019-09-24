.libPaths("/moto/ziab/users/yk2840/software/R_Libs/")

setwd("/moto/ziab/users/yk2840/Frog_workspace/HMM/Ancestry_nofilter_results/posterior/")

library(parallel)
numCores<-detectCores()

files<-list.files(pattern = "posterior$")


snp_keeps<-function(x, difffrac=.4){
  idx<-which.max(df[x,3:5])
  if(x==1){
    return(df[x,])}
  if(x!=1){
    if(idx!=which.max(df[x-1,3:5])){
      return(df[x,])
    }else{
    if(x %in% chrs_l){
      return(df[x,])
    }
    if(x==nrow(df)){
      return(df[x,])
    }else{
      if(df[x-1,1]==df[x,1]){
        if(abs(df[x-1,idx+2]-df[x,idx+2])>difffrac |
           abs(df[x+1,idx+2]-df[x,idx+2])>difffrac){
          return(df[x,])
          
        }}}}
  }
}

index<-NULL
for(fdf in files){
message(fdf)
df <- read.delim(fdf, stringsAsFactors=FALSE)
df$id<-1:nrow(df)


chrs_l<-order(df$chrom)[!duplicated(sort(df$chrom))] 

test<-do.call(rbind, mclapply(1:nrow(df), snp_keeps, mc.cores = numCores))
index<-unique(c(index, test$id))
}

write.table(index, paste0("snps_thinned_frac0.4.txt"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
