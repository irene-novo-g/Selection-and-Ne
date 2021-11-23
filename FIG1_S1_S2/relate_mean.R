
pop <- read.delim("popsizeall", header=FALSE)
colnames(pop) <- c("gen", "Ne")
pop <- subset(pop, (pop$gen >= 150) & (pop$gen <= 350) & (pop$Ne != Inf))

meangen <- vector()
for(i in 1:length(unique(pop$gen))){
  g <- pop$gen[i]
  gpop <- subset(pop, pop$gen == g)
  meangen[i] <- mean(gpop$Ne)
}

if(length(unique(pop$gen)) == length(meangen)){
  sink("NeRelate")
  cat("Mean_NeRelate_gens_150-350=\t")
  cat(mean(meangen))
  cat("\n")
  for(i in 1:length(meangen)){
    cat(pop$gen[i])
    cat("\t")
    cat(meangen[i])
    cat("\n")
  }
  sink()
} else { print("PRINTING ERROR: CHECK R SCRIPT")  }
