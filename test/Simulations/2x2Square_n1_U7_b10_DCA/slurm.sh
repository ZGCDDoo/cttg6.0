#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=2000
#SBATCH --account=def-tremblay


module reset
module load nixpkgs/16.09  gcc/5.4.0 armadillo boost-mpi

ITER=1
ITERMAX=10
myExe=cttg_DCA

if [ -a logfile ]
  then rm logfile
fi
rm tktilde.arma tloc.arma hybFM.arma config*.dat

while [ $ITER -le $ITERMAX ]
do
  echo begin iteration $ITER at: `date` >> logfile 

  srun $myExe params ${ITER}  

  echo end iteration $ITER at: `date` >> logfile
  ITER=$[$ITER+1]
done

