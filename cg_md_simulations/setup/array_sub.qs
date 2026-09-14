#!/bin/bash -l  

#SBATCH --nodes=1 
#SBATCH --ntasks=8
# SBATCH --partition=_workgroup_  
#SBATCH --partition=standard  
#SBATCH --job-name=charge.job #
#SBATCH --cpus-per-task=1  
#SBATCH --mem=2G  
#SBATCH --time=0-56:00:00
# SBATCH --time-min=0-00:10:00  #Do this
# SBATCH --output=slurms/%x-%j.out  
# SBATCH --error=slurms/%x-%j.out  
#SBATCH --mail-user='skronen@udel.edu'  
#SBATCH --mail-type=END,FAIL  
#SBATCH --export=NONE  
#UD_QUIET_JOB_SETUP=YES  
#UD_MACHINE_FILE_FORMAT='%h%[:]C'  
#export UD_JOB_EXIT_FN_SIGNALS="SIGTERM EXIT"  
#SBATCH --requeue
#SBATCH --array=1-N_ARRAY_TASKS_REPLACE


###INPUT PARAMETERS###
lmpfile="pgp_equil.in"
rs_file="./equil.rs.*"

FOLDER_FILE='./index.txt'
FOLDER=$(sed -n "$SLURM_ARRAY_TASK_ID p" "$FOLDER_FILE")
cd "$FOLDER"

###VALET### 
vpkg_require openmpi
. /opt/shared/slurm/templates/libexec/openmpi.sh  

if compgen -G "./post_equil.rs" > /dev/null; then
	echo "Skipping equil"
	
	else
        ###RUN### 
	mpirun lmp_2024_kspace -in $lmpfile
	mpi_rc=$?  

	if [ $mpi_rc -eq "0" ]; #if run finishes, 
	then rm $rs_file
	else exit 1
	fi
fi

###INPUT PARAMETERS###
lmpfile="pgp_prod.in"
rs_file="./prod.rs.*"

###DETERMINE IF RESTARTING### 
if compgen -G "./prod.rs.*" > /dev/null; then
	rsval=1 #if restart exists, set restart value to 1
else
	rsval=0
fi
 
###RUN### 
mpirun lmp_2024_kspace -in $lmpfile -v restarting $rsval
mpi_rc=$?  

if [ $mpi_rc -eq "0" ]; #if run finishes, 
then rm $rs_file
else exit 1
fi 

srun --cpus-per-task=1 bash filter_bonds.sh
mv filtered_pairs.txt pairs.txt


exit 0 

