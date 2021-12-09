#script_SLIM3
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

########################################################

module load gsl/2.1

#Working directory
WDIR=$PWD 

mkdir -p RESULTS_SLIM_$INPUT/

#Scratch directory
mkdir -p /state/partition1/IreneSlim3_$n/$SLURM_JOBID/

###################### TRANSFER TO SCRATCH ########################

cp $INPUT /state/partition1/IreneSlim3_$n/$SLURM_JOBID/
cp READSLIMOUT3_SLIM3 /state/partition1/IreneSlim3_$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/IreneSlim3_$n/$SLURM_JOBID

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

#cp -r /state/partition1/IreneSlim3_$n/$SLURM_JOBID/slimout $WDIR/RESULTS_SLIM_$INPUT/slimout$r

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

cp -r /state/partition1/IreneSlim3_$n/$SLURM_JOBID/timefile $WDIR/RESULTS_SLIM_$INPUT/
cp -r /state/partition1/IreneSlim3_$n/$SLURM_JOBID/chromosome$r.map $WDIR/RESULTS_SLIM_$INPUT/chromosome$r.map
cp -r /state/partition1/IreneSlim3_$n/$SLURM_JOBID/chromosome$r.ped $WDIR/RESULTS_SLIM_$INPUT/chromosome$r.ped

rm slimout
rm *.map
rm *.ped

######################## END OF REPLICATES #########################
done
######################## END OF REPLICATES #########################


######################## SCRATCH CLEANING #########################

rm $WDIR/$SLURM_JOBID.* 2> /dev/null
rm -rf /state/partition1/IreneSlim3_$n/$SLURM_JOBID/ 2> /dev/null
