#!/bin/sh
cd "emboss"
echo "[INFO] Running command: getorf -minsize $3 -sequence $1 -outseq $2"
getorf -minsize $3 -sequence $1 -outseq $2
