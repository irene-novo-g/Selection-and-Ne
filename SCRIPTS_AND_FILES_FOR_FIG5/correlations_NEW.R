
set.seed(3737)

require(dplyr)


#FINNISH GENERACIONES 25-50

load("Nereg2cM_FINNISH_GEN25-50_Nepmax100000.rda")
Nereg$Froh <- NULL
Nereg$Froh_mediumsize <- NULL
Nereg$F1 <- NULL
Nereg$F2 <- NULL
Nereg$F3 <- NULL
Nereg <- Nereg[c("chr", "reg", "Nef", "Nep", "heter", "pi", "nsnps", "npairs", "Fplink", "D", "Froh_big", "c", "skel", "nonskel", "dis", "gwas", "B", "nsnpsmaflow", "nsnpsmafhigh", "prop_maflow", "maf", "lof", "missense", "bp", "P", "skelbp", "nonskelbp", "disbp", "gwasbp", "lofbp", "missbp", "Ner", "Ner150", "Ner200", "genesbp")]

#KORYAKS GENERACIONES 25-50

load("C:/Users/user/OneDrive/TFM/Ne/Ne2cM_KORYAKS_hg19_GEN25-50/Nereg2cM_5c_fglobal_Nep100000.rda")
Nereg <- Nereg[c("chr", "reg", "Nef", "Nep", "het", "pi", "nsnps", "npairs", "Fplink", "D", "Froh_big", "rec", "skel", "nonskel", "dis", "gwas", "B", "nsnpsmaflow", "nsnpsmafhigh", "prop_maflow", "maf", "lof", "missense", "bp", "P", "skelbp", "nonskelbp", "disbp", "gwasbp", "lofbp", "missbp", "Ner", "Ner150", "Ner200", "genesbp")]

Pval = 0.05/2

metodo_cor = "spearman"


correlations.1variable_mod <- function(NUM.VARIABLE) {
  COR <- list()
  SIG <- list()
  P <- list()
  for(num.variable in 3:ncol(Nereg)){
    COR[[num.variable-2]] <- numeric(length = 23)
    SIG[[num.variable-2]] <- logical(length = 23)
    P[[num.variable-2]] <- numeric(length = 23)
    # 
    # for (chrom in c(1:22)) {
    #   chromosome <- dplyr::filter(Nereg, chr == chrom)
    #   if(sd(chromosome[[num.variable]], na.rm = T) == 0 || sd(chromosome[[NUM.VARIABLE]], na.rm = T) == 0){
    #     COR[[num.variable-2]][chrom] <- NA
    #     SIG[[num.variable-2]][chrom] <- NA
    #     P[[num.variable-2]][chrom] <- NA
    #     print("standard deviation is 0")
    #   } else {
    #     COR[[num.variable-2]][chrom] <- cor(chromosome[[NUM.VARIABLE]], chromosome[[num.variable]], method = metodo_cor, use = "na.or.complete")
    #     coresp <- vector()
    #     for (i in 1:10000) {
    #       creord <- sample(chromosome[[NUM.VARIABLE]], size = nrow(chromosome), replace = F)
    #       coresp[i] <- cor(chromosome[[num.variable]], creord, method = metodo_cor, use = "na.or.complete")
    #     }
    #     coresp <- sort(coresp, decreasing = T)
    #     maxcor <- coresp[1:(Pval*10000)]
    #     upperthreshold <- min(maxcor)
    #     mincor <- coresp[-(1:((1-Pval)*10000))]
    #     lowerthreshold <- max(mincor)
    #     
    #     
    #     if(COR[[num.variable-2]][chrom] > upperthreshold){
    #       SIG[[num.variable-2]][chrom] <- "HIGH"
    #     } else if (COR[[num.variable-2]][chrom] < lowerthreshold){
    #       SIG[[num.variable-2]][chrom] <- "LOW"
    #     } else {
    #       SIG[[num.variable-2]][chrom] <- NA
    #     }
    #     
    #     coresp <- sort(coresp)
    #     cordif <- abs((COR[[num.variable-2]][chrom])-coresp)
    #     P[[num.variable-2]][chrom] <- which.min(cordif)/10000
    #     
    #   }
    # }
    
    #GLOBAL
    chromosome <- dplyr::filter(Nereg, chr %in% 1:22)
    if(sd(chromosome[[num.variable]], na.rm = T) == 0 || sd(chromosome[[NUM.VARIABLE]], na.rm = T) == 0 || sum(!is.na(chromosome[[num.variable]])) == 0){
      COR[[num.variable-2]][23] <- NA
      SIG[[num.variable-2]][23] <- NA
      P[[num.variable-2]][23] <- NA
    } else {    
      COR[[num.variable-2]][23] <- cor(chromosome[[NUM.VARIABLE]], chromosome[[num.variable]], method = metodo_cor, use = "na.or.complete")
      corespgl <- vector()
      for (i in 1:10000) {
        creord <- sample(chromosome[[NUM.VARIABLE]], size = nrow(chromosome), replace = F)
        corespgl[i] <- cor(chromosome[[num.variable]], creord, method = metodo_cor, use = "na.or.complete")
      }
      corespgl <- sort(corespgl, decreasing = T)
      maxcor <- corespgl[1:(Pval*10000)]
      upperthreshold <- min(maxcor)
      mincor <- corespgl[-(1:((1-Pval)*10000))]
      lowerthreshold <- max(mincor)
    
      if(COR[[num.variable-2]][23] > upperthreshold){
        SIG[[num.variable-2]][23] <- "HIGH"
      } else if (COR[[num.variable-2]][23] < lowerthreshold){
        SIG[[num.variable-2]][23] <- "LOW"
      } else {
        SIG[[num.variable-2]][23] <- NA
      }
      
      corespgl <- sort(corespgl)
      cordifgl <- abs((COR[[num.variable-2]][23])-corespgl)
      P[[num.variable-2]][23] <- which.min(cordifgl)/10000
    }
    print(num.variable/ncol(Nereg))
  }
  
  correlations <- data.frame(Nef = COR[[1]], Nefsig = SIG[[1]], NefP = P[[1]],
                             Nep = COR[[2]], Nepsig = SIG[[2]], NepP = P[[2]], 
                             heter = COR[[3]], hetersig = SIG[[3]], heterP = P[[3]], 
                             pi = COR[[4]], pisig = SIG[[4]], piP = P[[4]],
                             nsnps = COR[[5]], nsnpssig = SIG[[5]], nsnpsP = P[[5]],
                             npairs = COR[[6]], npairssig = SIG[[6]], npairsP = P[[6]], 
                             Fplink = COR[[7]], Fplinksig = SIG[[7]], FplinkP = P[[7]],
                             D = COR[[8]], Dsig = SIG[[8]], DP = P[[8]],
                             Froh_big = COR[[9]], Froh_bigsig = SIG[[9]], Froh_bigP = P[[9]],
                             L = COR[[10]], Lsig = SIG[[10]], LP = P[[10]],
                             skel = COR[[11]], skelsig = SIG[[11]], skelP = P[[11]],
                             nonskel = COR[[12]], nonskelsig = SIG[[12]], nonskelP = P[[12]],
                             dis = COR[[13]], dissig = SIG[[13]], disP = P[[13]],
                             gwas = COR[[14]], gwassig = SIG[[14]], gwasP = P[[14]],
                             B = COR[[15]], Bsig = SIG[[15]], BP = P[[15]],
                             qlow = COR[[16]], qlowsig = SIG[[16]], qlowP = P[[16]],
                             qhigh = COR[[17]], qhighsig = SIG[[17]], qhighP = P[[17]],
                             prppqlow = COR[[18]], prppqlowsig = SIG[[18]], prppqlowP = P[[18]],
                             maf = COR[[19]], mafsig = SIG[[19]], mafP = P[[19]],
                             lof = COR[[20]], lofsig = SIG[[20]], lofP = P[[20]],
                             miss = COR[[21]], misssig = SIG[[21]], missP = P[[21]],
                             poly = COR[[23]], polysig = SIG[[23]], polyP = P[[23]],
                             skelbp = COR[[24]], skelbpsig = SIG[[24]], skelbpP = P[[24]],
                             nonskelbp = COR[[25]], nonskelbpsig = SIG[[25]], nonskelbpP = P[[25]],
                             disbp = COR[[26]], disbpsig = SIG[[26]], disbpP = P[[26]],
                             gwasbp = COR[[27]], gwasbpsig = SIG[[27]], gwasbpP = P[[27]],
                             lofbp = COR[[28]], lofbpsig = SIG[[28]], lofbpP = P[[28]],
                             missbp = COR[[29]], missbpsig = SIG[[29]], missbpP = P[[29]],
                             Nerel = COR[[30]], Nerelsig = SIG[[30]], NerelP = P[[30]],
                             Nerel150 = COR[[31]], Nerel150sig = SIG[[31]], Nerel150P = P[[31]],
                             Nerel200 = COR[[32]], Nerel200sig = SIG[[32]], Nerel200P = P[[32]],
                             genes = COR[[33]], genessig = SIG[[33]], genesP = P[[33]]
  )
  
  return(correlations)
}


#Nereg <- dplyr::filter(Nereg, chr != 6)

corNef <- correlations.1variable_mod(3) #Ne formula
write.csv(corNef, file = "corNef.csv")
corNep <- correlations.1variable_mod(4) #Ne program
write.csv(corNep, file = "corNep.csv")
#corh <- correlations.1variable_mod(5) #heterocigosis
#write.csv(corh, file = "corh.csv")
corpi <- correlations.1variable_mod(6) #pi
write.csv(corpi, file = "corpi.csv")
#cornsnps <- correlations.1variable_mod(7) #number of SNPs
#write.csv(cornsnps, file = "cornsnps.csv")
#cornpairs <- correlations.1variable_mod(8) #number of pairs
#write.csv(cornpairs, file = "cornpairs.csv")
corFplink <- correlations.1variable_mod(9) #F plink
write.csv(corFplink, file = "corFplink.csv")
corD <- correlations.1variable_mod(10) #D
write.csv(corD, file = "corD.csv")
corFroh_big <- correlations.1variable_mod(11) #F roh big (rohs >= 100Kb)
write.csv(corFroh_big, file = "corFroh_big.csv")
corL <- correlations.1variable_mod(12) #L (recombination rate in cM/Mb)
write.csv(corL, file = "corL.csv")
#corskel <- correlations.1variable_mod(13) #skeletal SNPs
#write.csv(corskel, file = "corskel.csv")
#cornonskel <- correlations.1variable_mod(14) #non-skeletal SNPs
#write.csv(cornonskel, file = "cornonskel.csv")
#cordis <- correlations.1variable_mod(15) #disease SNPs
#write.csv(cordis, file = "cordis.csv")
#corgwas <- correlations.1variable_mod(16) #gwas SNPs
#write.csv(corgwas, file = "corgwas.csv")
corB <- correlations.1variable_mod(17) #B
write.csv(corB, file = "corB.csv")
#corqlow <- correlations.1variable_mod(18) #number of SNPs whose MAF <=0.05
#write.csv(corqlow, file = "corqlow.csv")
#corqhigh <- correlations.1variable_mod(19) #number of SNPs whose MAF > 0.05
#write.csv(corqhigh, file = "corqhigh.csv")
#corpropqlow <- correlations.1variable_mod(20) #(number of SNPs whose MAF <=0.05) / number of SNPs
#write.csv(corpropqlow, file = "corpropqlow.csv")
cormaf <- correlations.1variable_mod(21) # mean MAF
write.csv(cormaf, file = "cormaf.csv")
#corlof <- correlations.1variable_mod(22) #number of LOF variants according to ExAC.r0.3.1.sites.vep.vcf
#write.csv(corlof, file = "corlof.csv")
#cormiss <- correlations.1variable_mod(23) #number of missense variants according to ExAC.r0.3.1.sites.vep.vcf
#write.csv(cormiss, file = "cormiss.csv")
#base pairs
corpoly <- correlations.1variable_mod(25) #polymorphism: nºSNPs/pb
write.csv(corpoly, file = "corpoly.csv")
#corskelbp <- correlations.1variable_mod(26) #skel/pb
#write.csv(corskelbp, file = "corskelbp.csv")
#cornonskelbp <- correlations.1variable_mod(27) #nonskel/pb
#write.csv(cornonskelbp, file = "cornonskelbp.csv")
#cordisbp <- correlations.1variable_mod(28) #dis/pb
#write.csv(cordisbp, file = "cordisbp.csv")
#corgwasbp <- correlations.1variable_mod(29) #gwas/pb
#write.csv(corgwasbp, file = "corgwasbp.csv")
corlofbp <- correlations.1variable_mod(30) #lof/pb
write.csv(corlofbp, file = "corlofbp.csv")
cormissbp <- correlations.1variable_mod(31) #missense/pb
write.csv(cormissbp, file = "cormissbp.csv")
#corner <- correlations.1variable_mod(32) #Ne Relate recent generations
#write.csv(corner, file = "corner.csv")
#corner150 <- correlations.1variable_mod(33) #Ne Relate gens 150-200
#write.csv(corner150, file = "corner150.csv")
#corner200 <- correlations.1variable_mod(34) #Ne Relate gens 200-250
#write.csv(corner200, file = "corner200.csv")
corgenedens <- correlations.1variable_mod(35) #genedensity
write.csv(corgenedens, file = "corgenedens.csv")

