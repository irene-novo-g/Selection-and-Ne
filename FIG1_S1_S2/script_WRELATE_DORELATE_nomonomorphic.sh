#script_WRELATE_DORELATE_nomonomorphic
#!/bin/bash
#$ -cwd

#rm script_WRELATE_DORELATE_nomonomorphic.sh.*

################## INPUT FILE AND REPLICATES AS ARGUMENTS ############

module load gcc/7.2.0
module load gsl/2.1
module load R/3.3.0
module load curl/7.47.1
module load glibc/2.14


################## INPUT FILE AND REPLICATES AS ARGUMENTS ############

#Check number of arguments
if [ $# -ne 3 ]  
then
	echo "Usage: $0 <NIND> <NREPS> <n>" 
	exit 1
fi

#Set arguments
NIND=$1
REPS=$2
n=$3

########################################################

#Working directory
WDIR=$PWD 

mkdir -p OUT_RELATE/

#Scratch directory
mkdir -p /state/partition1/IreneRelate$n/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

###################### TRANSFER TO SCRATCH ########################

cp READMAPPED_WRITERELATE_SIMS /state/partition1/IreneRelate$n/$SLURM_JOBID/
cp relate_mean.R /state/partition1/IreneRelate$n/$SLURM_JOBID/
cp relate_delete_monomorphic.R /state/partition1/IreneRelate$n/$SLURM_JOBID/
cp -r /home/invitado1/IRENE/RELATE /state/partition1/IreneRelate$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/IreneRelate$n/$SLURM_JOBID

################################ REPLICATES #############################

for ((chr=1; chr<=REPS; chr++))
do

cp $WDIR/chromosome$chr.map /state/partition1/IreneRelate$n/$SLURM_JOBID/
cp $WDIR/chromosome$chr.ped /state/partition1/IreneRelate$n/$SLURM_JOBID/

###################### READMAPPED_WRITERELATE #########################

START=$(date +%s)

if [ -s "chromosome$chr.map" ]; then

awk '{print $1,$3,$4}' chromosome$chr.map > dataBP.map  # remove name of SNPs ($2), hard to read in READMAPPED_WRITERELATE and irrelevant
awk '{$1=$2=$3=$4=""; print $0}' chromosome$chr.ped > dataBP.ped # remove name of population, of individual, and 2 first numbers. remain sex data

./READMAPPED_WRITERELATE_SIMS<<@
$NIND	N
$chr	NCHR
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "READMAPPED_WRITERELATE took 		$DIFF seconds" >> timefile
echo "chr = $chr" >> timefile

rm dataBP.map
rm dataBP.ped
cp data.haps chromosome_relate$chr.haps
cp data.sample chromosome_relate$chr.sample
cp data.map chromosome_relate$chr.map

cp -r /state/partition1/IreneRelate$n/$SLURM_JOBID/timefile $WDIR/OUT_RELATE/

#cp -r /state/partition1/IreneRelate$n/$SLURM_JOBID/chromosome_relate$chr.haps $WDIR/OUT_RELATE/
#cp -r /state/partition1/IreneRelate$n/$SLURM_JOBID/chromosome_relate$chr.map $WDIR/OUT_RELATE/
#cp -r /state/partition1/IreneRelate$n/$SLURM_JOBID/chromosome_relate$chr.sample $WDIR/OUT_RELATE/
#cp -r /state/partition1/IreneRelate$n/$SLURM_JOBID/data.poplabels $WDIR/OUT_RELATE/

Rscript --vanilla relate_delete_monomorphic.R $chr

#cp -r /state/partition1/IreneRelate$n/$SLURM_JOBID/chromosome_relate_nomonomorphic$chr.haps $WDIR/OUT_RELATE/


#################### RELATE ###########################

START=$(date +%s)

HAPS=$((NIND*2)) 

/state/partition1/IreneRelate$n/$SLURM_JOBID/RELATE/bin/Relate \
      --mode All \
      -m 1.00e-8 \
      -N $HAPS \
      --haps chromosome_relate_nomonomorphic$chr.haps \
      --sample chromosome_relate$chr.sample \
      --map chromosome_relate$chr.map \
      --seed 1 \
      -o out$chr

/state/partition1/IreneRelate$n/$SLURM_JOBID/RELATE/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i out$chr \
              -m 1.00e-8 \
              --poplabels data.poplabels \
	      --years_per_gen 80 \
	      --num_bins 300 \
	      --threshold 50 \
              --seed 1 \
              -o out_popsize$chr

rm *.anc
rm *.mut
rm *.

#cp -r /state/partition1/IreneRelate$n/$SLURM_JOBID/out_popsize$chr.coal $WDIR/OUT_RELATE/

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "RELATE took 		$DIFF seconds" >> timefile
echo "chr = $chr" >> timefile


##################### GENERATE READABLE RESULTS ####################################################################

if [ -s "out_popsize$chr.coal" ]; then
	echo "out_popsize$chr.coal exists"
    sed '1d' out_popsize$chr.coal > out_popsize${chr}_m.coal

	awk '
	{ 
	    for (i=1; i<=NF; i++)  {
	        a[NR,i] = $i
 	   }
	}
	NF>p { p = NF }
	END {    
 	   for(j=1; j<=p; j++) {
  	      str=a[1,j]
 	       for(i=2; i<=NR; i++){
	            str=str" "a[i,j];
 	       }
 	       print str
 	   }
	}' out_popsize${chr}_m.coal > out_popsize${chr}_trans.coal

	awk '{ print $1,"\t", ($2^-1)*0.5}' out_popsize${chr}_trans.coal > out_popsize${chr}_NE.coal


cp out_popsize${chr}_NE.coal popsize$chr

cat popsize$chr >> popsizeall

Rscript --vanilla relate_mean.R

#cp /state/partition1/IreneRelate$n/$SLURM_JOBID/popsize$chr $WDIR/OUT_RELATE/
cp -r /state/partition1/IreneRelate$n/$SLURM_JOBID/popsizeall $WDIR/OUT_RELATE/
cp /state/partition1/IreneRelate$n/$SLURM_JOBID/NeRelate $WDIR/OUT_RELATE/

fi

rm popsize$chr
rm *.coal
rm *.rate
rm *.bin
rm *.anc
rm *.dist
rm *.mut
rm *.pdf
rm chromosome$chr.map
rm chromosome$chr.ped
rm chromosome_relate_nomonomorphic$chr.haps
rm chromosome_relate$chr.haps
rm chromosome_relate$chr.sample
rm chromosome_relate$chr.map


fi

done
######################## END OF REPLICATES #########################

######################## SCRATCH CLEANING #########################

rm $WDIR/$SLURM_JOBID.* 2> /dev/null
rm -rf /state/partition1/IreneRelate$n/$SLURM_JOBID/ 2> /dev/null
