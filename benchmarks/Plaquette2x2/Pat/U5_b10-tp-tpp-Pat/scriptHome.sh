
SC_DIR=../
IS_DIR=../../../impuritySolver/


ITER=1
ITERMAX=100

if [ -a logfile ]
  then rm logfile
fi

if [ $ITER -eq 0 ]
then
  $SC_DIR/CDMFT params ${ITER}

  ITER=$[$ITER+1]
fi

while [ $ITER -le $ITERMAX ]
do
  

  echo begin iteration $ITER at: `date` >> logfile 
  $SC_DIR/GA params${ITER}
  
  echo start impurity solver at: `date` >> logfile
  $IS_DIR/IS params${ITER}
  echo end impurity solver at: `date` >> logfile
  
  echo begin self-consistency at: `date` >> logfile
  $SC_DIR/CDMFT params ${ITER}
  echo end self-consistency at: `date` >> logfile
  
  echo end iteration $ITER at: `date` >> logfile
  ITER=$[$ITER+1]
done

