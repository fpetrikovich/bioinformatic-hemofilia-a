#!/bin/sh
cd $1
echo "[INFO] RUNNING COMMAND: getorf -minsize $5 -sequence $2 -outseq $3"
getorf -minsize $5 -sequence $2 -outseq $3
echo "[INFO] RUNNING COMMAND: patmatmotifs -sequence $3 -outfile $4 -full"
patmatmotifs -sequence $3 -outfile $4 -full