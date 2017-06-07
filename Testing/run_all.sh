#!/usr/bin/sh

N=50
MAX_ITER=5000
EPS=0.01
nw=256

RES_FF_FILE="res_ff.csv"
RES_SEQ_FILE="res_seq.csv"

echo "Running sequential"
echo "nw, it, time"  > $RES_SEQ_FILE
JacobiSeq/cmake-build-debug/JacobiSeq $N $MAX_ITER $EPS >> $RES_SEQ_FILE

echo "Running fast-flow"
echo "nw, it, time"  > $RES_FF_FILE
for ((i = 2; i <= $nw; i+=5))
do
    echo "Working with $i"
    JacobiFF/cmake-build-debug/JacobiFF $N $MAX_ITER $EPS $i >> $RES_FF_FILE
done
