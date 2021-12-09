
pop <- read.delim("popsizeall", header=FALSE)
colnames(pop) <- c("gen", "Ne")
pop <- subset(pop, pop$Ne != Inf)

meangen <- vector()
for(i in 1:length(unique(pop$gen))){
  g <- sort(unique(pop$gen))[i]
  gpop <- subset(pop, pop$gen == g)
  meangen[i] <- mean(gpop$Ne)
  print(g)
  print(meangen[i])
}

if(length(unique(pop$gen)) == length(meangen)){
  sink("NeRelate")
  for(i in 1:length(meangen)){
    cat(sort(unique(pop$gen))[i])
    cat("\t")
    cat(meangen[i])
    cat("\n")
  }
  sink()
} else { print("PRINTING ERROR: CHECK R SCRIPT")  }