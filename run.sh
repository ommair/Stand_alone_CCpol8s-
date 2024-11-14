#!/bin/bash
#
#SBATCH --job-name=ccpol23
#SBATCH --output=myjob-%A_%a.out
#SBATCH --error=myjob-%A_%a.err
#SBATCH --partition=kolos2
#
###SBATCH --constraint="infiniband"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=14-00:00:0
#SBATCH --mem-per-cpu=1024MB

echo "***********************************************"
hostname
ulimit -a

##SBATCH --ntasks-per-core=1
##SBATCH --exclusive


# this job requests 32 cores. Cores can be selected from various nodes.

df -h 
ls 
echo $JOBDIR

module load openmpi/3.0.0-gcc 

time /data/ommair/My_MD_MC_Codes/MD_FCC_010222_VV_NBind_ES_3B_restart_old/execute

#EXE=/data/ommair/My_MD_MC_Codes/MD_FCC_030220_VV_changes/execute

#WRKDIR=/scratch/local/ommair/$SLURM_JOB_ID
#ILEDIR=$PWD


#export OMP_NUM_THREADS=4

#mkdir -p $WRKDIR


#cp $ILEDIR/* $WRKDIR/
#cp $ILEDIR/STATIS $WRKDIR/

#cd $WRKDIR/
#echo $HOSTNAME

#hostname
#which mpirun 

#mpirun -np ${OMP_NUM_THREADS} $EXE  
#mpirun -np 1 $EXE
#rm slurm*.out
#cp config_dlpoly.xyz $ILEDIR/
#cp OUTPUT $ILEDIR/
#cp HISTORY.xyz $ILEDIR/

