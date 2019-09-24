.libPaths("/moto/ziab/users/yk2840/software/R_Libs")

library(reshape2)
library(stringr)

#setwd("/moto/ziab/users/yk2840/Frog_workspace/HMM/Ancestry_nofilter_results/")
setwd("/moto/ziab/users/yk2840/Frog_workspace/HMM/Ancestry_nofilter_results/raw_allelecounts/")
load("/moto/ziab/users/yk2840/Frog_workspace/scripts_ahmm/ancestry_panel.RData")
load("/moto/ziab/users/yk2840/Frog_workspace/scripts_ahmm/fixed_pet_variants_allelecounter.RData")

message("loading")
files<-list.files(pattern = "_ac_nofilter.txt")

library(parallel)
numcores<-detectCores()

ac_combine<-function(x){
  
  message(x)
  samp <- read.delim(x, 
                     header=FALSE, comment.char="#", 
                     stringsAsFactors=FALSE)
  colnames(samp)<-c("CHROM", "POS", "A", "C","G", "T", "Depth")
  
  
  samp<-melt(samp, id.vars = colnames(samp)[c(1,2,7)])
  samp$ID<-gsub(".txt", "", basename(x))
  samp$UID<-with(samp, paste0(CHROM, ":", POS, "_", variable))
  
  df<-data.frame(REF = samp$value[match(fixed_pet_variants_allelecounter$UID,
                                        samp$UID)],
                 ALT = samp$value[match(fixed_pet_variants_allelecounter$UAID,
                                        samp$UID)])
  colnames(df)<-c(paste0(str_split_fixed(x, "_", 2)[,1], "_R"),
                  paste0(str_split_fixed(x, "_", 2)[,1], "_A"))
  return(df) 
}

totals <- mclapply(files, ac_combine, mc.cores = numcores)
head(totals)
totals <- do.call(cbind, totals)
head(totals)

# Reduce analysis by removing all zeros:
remove.idx<-apply(totals, 1, sum)
remove.idx<-which(remove.idx==0)


write.table(totals, file = paste0(gsub("-", "_", Sys.Date()), "_sample_panel_notfiltered.txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

totals<-totals[remove.idx,]
write.table(totals, file = paste0(gsub("-", "_", Sys.Date()), "_sample_panel_filtered_norecomb.txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

ancestry_panel<-ancestry_panel[remove.idx,]
#ancestry_panel$cM<-c(0,diff(ancestry_panel$V2)*10^-8)/1000
ancestry_panel$cM<-c(0,diff(ancestry_panel$V2))


totals<-cbind(ancestry_panel, totals)
write.table(totals, file = paste0(gsub("-", "_", Sys.Date()), "_ac_panel_norecomb.txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)

library(stringr)
ploidy<-data.frame(name = str_split_fixed(files, "_", 2)[,1],
                   ploidy = 2)
write.table(ploidy, file = paste0(gsub("-", "_", Sys.Date()), "_ac_ploidy.txt"), 
            quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)


               
