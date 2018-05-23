

#!/bin/bash
iter=1
iterMax=10
myExe=cdmft_triangle2x2_aux
nprocess=4

rm tktilde.arma tloc.arma hybFM.arma config.dat
while [  $iter -lt $iterMax ];
    do
        mpirun -np $nprocess $myExe params $iter
	    iter=$[iter+1]
    done


