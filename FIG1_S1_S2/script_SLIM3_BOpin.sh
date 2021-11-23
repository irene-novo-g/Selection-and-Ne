#script_SLIM3_BOpin
#!/bin/bash
#$ -cwd

#rm script_SLIM3_BOpin.sh.*

################## INPUT FILE AND REPLICATES AS ARGUMENTS ############

#Check number of arguments
if [ $# -ne 5 ]  
then
	echo "Usage: $0 <slimINPUTfile> <NIND> <cMMb> <REPS> <n>" 
	exit 1
fi

#Set arguments
INPUT=$1
NIND=$2
cMMb=$3
REPS=$4
n=$5

### Other parameters

PHASE=2
DIST=0
NGEN=1000 ### Number of generations for which linkage data is obtained in bins
NBIN=200  ### Number of bins (e.g. 1000, so that each bin includes NGEN/NBIN generations)
SNP=-99   ### Maximum number of SNPs to be used per chromosome (if SNP=-99 or PHASE=1 it is not applied)
MAF=0.0   ### Minimum allele frequency (0-1)
ZERO=1    ### 0: Remove SNPs with zeroes (1: allow for them)
SEPARATE=1

########################################################

module load gsl/2.1

#Working directory
WDIR=$PWD 

mkdir -p RESULTS_SLIM_BOpin$INPUT/

#Scratch directory
mkdir -p /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

###################### TRANSFER TO SCRATCH ########################

cp $INPUT /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/
cp READSLIMOUT3_SLIM3 /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/
cp LD_SNP_REAL_BOpin3 /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/
cp BOpin /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/
cp SUMM_REP_BOpin /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID

################################ REPLICATES #############################

module load SLiM/3.3.2

for ((r=1; r<=$REPS; r++))
do

################################ SLIM #############################

START=$(date +%s)
slim $INPUT > slimout
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Slim3 took 		$DIFF seconds" >> timefile

cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/slimout $WDIR/RESULTS_SLIM_BOpin$INPUT/slimout$r

###################### READSLIMOUT2 #########################

START=$(date +%s)

./READSLIMOUT3_SLIM3<<@
$NIND	N
100000000  genomelen
1	NCHR
$cMMb	cMMb
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "READSLIMOUT3_SLIM3 took 		$DIFF seconds" >> timefile
echo "rep = $r" >> timefile

cp dataBP.map chromosome$r.map
cp dataBP.ped chromosome$r.ped

cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/timefile $WDIR/RESULTS_SLIM_BOpin$INPUT/
cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/chromosome$r.map $WDIR/RESULTS_SLIM_BOpin$INPUT/chromosome$r.map
cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/chromosome$r.ped $WDIR/RESULTS_SLIM_BOpin$INPUT/chromosome$r.ped

###### LD_SNP_REAL_BOpin #######

START=$(date +%s)

mv chromosome$r.map data.map
mv chromosome$r.ped data.ped

awk '{print $3 " " $4}' data.map > kk
mv kk data.map

./LD_SNP_REAL_BOpin3<<@
$NIND	NIND
$MAF		MAF
$PHASE	PHASE
$NGEN	NGEN 
$NBIN	NBIN
$ZERO	ZERO
$DIST	DIST
$cMMb	cMMb
$SEPARATE	SEPARATE
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Replicate $r took 		$DIFF seconds" >> timefile

mv outfileLD outfileLD.$r
mv parameters parameters.$r
mv heterozygosity heterozygosity.$r
cat heterozygosity.$r >> HETEROZYGOSITY_PI

cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/timefile $WDIR/RESULTS_SLIM_BOpin$INPUT/
##cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/outfileLD.$r $WDIR/RESULTS_SLIM_BOpin$INPUT/
##cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/parameters.$r $WDIR/RESULTS_SLIM_BOpin$INPUT/
##cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/heterozygosity.$r $WDIR/RESULTS_SLIM_BOpin$INPUT/
cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/HETEROZYGOSITY_PI $WDIR/RESULTS_SLIM_BOpin$INPUT/

######################## SUMM_REP_CHROM #########################

rm CHROM

for ((k=1; k<=r; k++))
do
cat outfileLD.$k >> CHROM
done

./SUMM_REP_BOpin<<@
$NGEN	NGEN
$NBIN	NBIN
$r		REPS
@

cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/outfileLD $WDIR/RESULTS_SLIM_BOpin$INPUT/

cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/list_allsnps $WDIR/RESULTS_SLIM_BOpin$INPUT/

############################# BOpin.cpp ##########################

./BOpin outfileLD

cp -r /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/outfileLD.BOpinNe $WDIR/RESULTS_SLIM_BOpin$INPUT/


######################## END OF REPLICATES #########################
done
######################## END OF REPLICATES #########################


######################## SCRATCH CLEANING #########################

rm $WDIR/$SLURM_JOBID.* 2> /dev/null
rm -rf /state/partition1/IreneSlimRepBOPIN$n/$SLURM_JOBID/ 2> /dev/null
