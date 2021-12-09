#!/bin/bash

module load cesga/2018
module load gcc/6.4.0
module load gsl/2.5
module load R/3.5.3
module load curl/7.61.1
module load glibc/2.28

cd $PWD

for((rep=1;rep<=100;rep+=10))
do

./msmc2_linux64bit --fixedRecombination -t 1 -o out$rep chromosome_msmc$rep chromosome_msmc$((rep+1)) chromosome_msmc$((rep+2)) chromosome_msmc$((rep+3)) chromosome_msmc$((rep+4)) chromosome_msmc$((rep+5)) chromosome_msmc$((rep+6)) chromosome_msmc$((rep+7)) chromosome_msmc$((rep+8)) chromosome_msmc$((rep+9))

awk '{ print $2 / 1e-8,"\t",  $3 / 1e-8, "\t", ($4 ^(-1))/(2e-8)}' out$rep.final.txt > out$rep.NE
awk '{ print ($1+$2)/2"\t",  $3}' out$rep.NE > out$rep.gNE

cut out$rep.gNE -f 1 > outgen$rep
paste outgen* > outgen

cut out$rep.gNE -f 2 > outne$rep
paste outne* > outne

Rscript --vanilla msmc_write_meangenNe.R

paste meangen meanne > globalNE

rm -r *.log
rm -r *.loop

done


