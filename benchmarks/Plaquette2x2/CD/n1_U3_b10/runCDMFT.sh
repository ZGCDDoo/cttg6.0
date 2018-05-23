

#!/bin/bash
iter=1
iterMax=10
myExe=cdmft_square2x2
nprocess=4

rm *Rank* tktilde.arma tloc.arma hybFM.arma config.dat
while [  $iter -lt $iterMax ];
    do
        mpirun -np $nprocess $myExe params $iter
	    iter=$[iter+1]
    done


