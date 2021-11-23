outgen <- read.delim("outgen", header=FALSE)
outgen[,1] <- NULL
meangen <- vector()
for(i in 1:nrow(outgen)){
  meangen[i] <- mean(as.numeric(outgen[i,]))
}
sink("meangen")
for(i in 1:length(meangen)){
  cat(meangen[i])
  cat("\n")
}
sink()



outne <- read.delim("outne", header=FALSE)
outne[,1] <- NULL
meanne <- vector()
for(i in 1:nrow(outne)){
  meanne[i] <- mean(as.numeric(outne[i,]))
}
sink("meanne")
for(i in 1:length(meanne)){
  cat(meanne[i])
  cat("\n")
}
sink()
