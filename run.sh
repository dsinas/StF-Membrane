#!/bin/bash

L=100
J=0.50
g=1.00

mkdir "data"
cd data

#for L in 200 180 160 140 120 100
#do
  LC_NUMERIC=en_US.UTF-8
  gseq=$(seq 0.98 0.02 1.02)
  for g in $gseq
  do
    dir="L"$L"-J"$J"-g"$g
    mkdir $dir
    cd $dir
    ../../stfMem -L $L -J $J -g $g
    cd ..
  done
#done

