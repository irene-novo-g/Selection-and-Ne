#script_WMSMC
#!/bin/bash
#$ -cwd

#rm script_WMSMC.sh.*

################## INPUT FILE AND REPLICATES AS ARGUMENTS ############

#Check number of arguments
if [ $# -ne 2 ]  
then
	echo "Usage: $0 <REPS> <n>" 
	exit 1
fi

#Set arguments
REPS=$1
n=$2

########################################################

module load gsl/2.1

#Working directory
WDIR=$PWD 

mkdir -p OUT_MSMC/

#Scratch directory
mkdir -p /state/partition1/IreneMSMC$n/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

###################### TRANSFER TO SCRATCH ########################

cp READMAPPED_WRITEMSMC /state/partition1/IreneMSMC$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/IreneMSMC$n/$SLURM_JOBID

################################ REPLICATES #############################

for ((chr=1; chr<=$REPS; chr++))
do

cp $WDIR/chromosome$chr.map /state/partition1/IreneMSMC$n/$SLURM_JOBID/
cp $WDIR/chromosome$chr.ped /state/partition1/IreneMSMC$n/$SLURM_JOBID/

###################### READMAPPED_WRITEMSMC #########################

if [ -s "chromosome$chr.map" ]; then

awk '{print $1,$3,$4}' chromosome$chr.map > dataBP.map #remove name of SNPs ($2), hard to read
awk '{$1=$2=$3=$4=""; print $0}' chromosome$chr.ped > dataBP.ped #remove name of population, indiviual and first 2 numbers. remain sex data

START=$(date +%s)

./READMAPPED_WRITEMSMC<<@
$chr	NCHR
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "READMAPPED_WRITEMSMC took 		$DIFF seconds" >> timefile
echo "rep = $chr" >> timefile

cp datamsmc chromosome_msmc$chr

cp -r /state/partition1/IreneMSMC$n/$SLURM_JOBID/timefile $WDIR/OUT_MSMC/
cp -r /state/partition1/IreneMSMC$n/$SLURM_JOBID/chromosome_msmc$chr $WDIR/OUT_MSMC/

rm chromosome_msmc$chr
rm chromosome$chr.map
rm chromosome$chr.ped

fi

######################## END OF REPLICATES #########################
done

######################## SCRATCH CLEANING #########################

rm $WDIR/$SLURM_JOBID.* 2> /dev/null
rm -rf /state/partition1/IreneMSMC$n/$SLURM_JOBID/ 2> /dev/null
