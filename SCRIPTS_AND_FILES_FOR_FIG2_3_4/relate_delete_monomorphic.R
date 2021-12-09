args=commandArgs(TRUE)
chr=args[1]

haps <- read.table(paste("chromosome_relate", chr, ".haps", sep = ""), quote="\"", comment.char="")
  
rows.to.delete <- vector()
for (i in 1:nrow(haps)) {
  row <- haps[i, 6:ncol(haps)]
  row <- as.numeric(row)
  rep <- length(unique(row))
  if (rep == 1) {  #MONOMORPHIC SNP
    rows.to.delete <- c(rows.to.delete, i)
  }
}

if(length(rows.to.delete) > 0){
  haps <- haps[-rows.to.delete,]
}

haps[,5] <- "T"
haps[,2] <- as.character(haps[,2])
haps[,4] <- "A"
  
sink(paste("chromosome_relate_nomonomorphic", chr, ".haps", sep = ""))
for(i in 1:nrow(haps)) {
  for(j in 1:ncol(haps)){
    cat(haps[i,j])
    cat(" ")
  }
  cat("\n")
}
sink()

