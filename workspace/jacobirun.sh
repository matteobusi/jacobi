#!/bin/bash

# JACOBI SETTINGS
MAX_ITER=100
EPS=0.001
N_VALS=(5000 10000 15000 30000)

# WORKERS SETTINGS
if [[ "$1" == "MIC" ]]
then
    echo -n "Running on MIC"
    nw=240
    WORKER_STEPS=16
else
    echo -n "Running on HOST"
    nw=16
    WORKER_STEPS=4
fi

echo " with $nw max. workers"

for N in "${N_VALS[@]}"
do    
    # FILES
    if [[ "$1" == "MIC" ]]
    then
	RES_FF_FILE="results/es_ff_mic_$N.csv"
	RES_SEQ_FILE="results/res_seq_mic_$N.csv"
	RES_PT_FILE="results/res_pt_mic_$N.csv"
    else
	RES_FF_FILE="results/res_ff_host_$N.csv"
	RES_SEQ_FILE="results/res_seq_host_$N.csv"
	RES_PT_FILE="results/res_pt_host_$N.csv"
    fi
    echo "N = $N"
    
    echo -n "Running sequential ... "
    echo "nw, it, time, err"  > $RES_SEQ_FILE
    if [[ "$1" == "MIC" ]]
    then
	ssh mic1 "bin/jacobim $N $MAX_ITER $EPS sq" >> $RES_SEQ_FILE
    else
	bin/jacobix $N $MAX_ITER $EPS sq >> $RES_SEQ_FILE
    fi
    echo "Done"
    
    echo "Running fast-flow"
    echo "nw, it, time, err"  > $RES_FF_FILE
    for i in 1 2 $(seq 4 $WORKER_STEPS $nw)
    do
	echo -n "Working with $i ... "
	if [[ "$1" == "MIC" ]]
	then
	    ssh mic1 "bin/jacobim $N $MAX_ITER $EPS ff $i" >> $RES_FF_FILE
	else
	    bin/jacobix $N $MAX_ITER $EPS ff $i >> $RES_FF_FILE
	fi
	echo "Done"
    done
	     
    echo "Running pthread"
    echo "nw, it, time, err" > $RES_PT_FILE
    for i in 1 2 $(seq 4 $WORKER_STEPS $nw)
    do
	echo -n "Working with $i ... "
        if [[ "$1" == "MIC" ]]
	then
	    ssh mic1 "bin/jacobim $N $MAX_ITER $EPS pt $i" >> $RES_PT_FILE
	else
	    bin/jacobix $N $MAX_ITER $EPS pt $i >> $RES_PT_FILE
        fi
	echo "Done"
    done
done
